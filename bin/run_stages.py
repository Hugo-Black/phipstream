#!/usr/bin/env python3
"""Run the pipeline stages in order from a single configuration file.

The stages are otherwise separate programs. On a batch scheduler they need to be
one command, because a job runs one command and the parameters have to be
recorded somewhere the job can read.

The configuration is an INI file with one [phipstream] section. Keys match the
stage options with dashes replaced by underscores. Trimming is skipped when no
adapters key is given, which suits a sample table that already points at trimmed
reads.

An optional [workflow] section holds parameters for the bundled workflow. With
--workflow those are rendered into a Nextflow config, taking the input tables and
adapters from [phipstream] so the two routes cannot disagree about a dataset.

With --submit the stages are not run. A job script is written that runs this
same command on a compute node, and it is submitted when qsub is on PATH.
"""
import argparse
import configparser
import csv
import os
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

BIN = Path(__file__).resolve().parent

# Stage order, the script each one runs, and the output that marks it complete.
STAGES = ("fastqc", "trim", "align", "score", "prioritise")
SCRIPTS = {"fastqc": "run_fastqc.py", "trim": "trim_reads.py",
           "align": "align_reads.py", "score": "call_enrichment.py",
           "prioritise": "prioritise_peptides.py"}
# Score matrix and hit matrix each scoring method leaves behind.
SCORE_FILES = {"beer": ("beer_posterior.csv.gz", "beer_hits.csv.gz"),
               "edger": ("edger_logpval.csv.gz", "edger_hits.csv.gz"),
               "larman_gp": ("larman_gp_mlxp.csv.gz", "larman_gp_hits.csv.gz"),
               "xu_zigp": ("xu_zigp_mlxp.csv.gz", "xu_zigp_hits.csv.gz")}


def load_config(path):
    parser = configparser.ConfigParser()
    # Keys are passed through as workflow parameter names, which are case
    # sensitive. Without this, run_BEER would arrive as run_beer and be ignored.
    parser.optionxform = str
    read = parser.read(path)
    if not read:
        sys.exit(f"cannot read config: {path}")
    if not parser.has_section("phipstream"):
        sys.exit(f"{path} needs a [phipstream] section")

    def section(name):
        return {k: v.strip() for k, v in parser.items(name) if v.strip()}

    workflow = section("workflow") if parser.has_section("workflow") else {}
    return section("phipstream"), workflow


def enabled(cfg, key):
    """True when a config value reads as a switch that is on."""
    return cfg.get(key, "").lower() in ("true", "yes", "1")


def require(cfg, key, path):
    if key not in cfg:
        sys.exit(f"{path} is missing the {key} key")
    return cfg[key]


def build_commands(cfg, config_path, out_dir):
    """Return [(stage, argv, output that marks it complete)] in run order."""
    sample_table = require(cfg, "sample_table", config_path)
    peptide_table = require(cfg, "peptide_table", config_path)
    mode = cfg.get("mode", "se_r1")
    preset = cfg.get("preset", "local")
    method = cfg.get("method", "beer")
    if method not in SCORE_FILES:
        sys.exit(f"method must be one of {', '.join(sorted(SCORE_FILES))}")

    trimmed_table = out_dir / "sample_table_trimmed.csv"
    counts = out_dir / "alignment" / "counts" / f"{mode}.{preset}.csv"
    score_name, hits_name = SCORE_FILES[method]
    posterior = out_dir / "enrichment" / score_name
    hits = out_dir / "enrichment" / hits_name

    plan = []
    if enabled(cfg, "fastqc"):
        fastqc_dir = out_dir / "fastqc"
        plan.append(("fastqc", [
            "--sample-table", sample_table,
            "--out-dir", str(fastqc_dir),
            "--threads", cfg.get("jobs", "4"),
        ], fastqc_dir / "fastqc_files.csv"))

    if "adapters" in cfg:
        plan.append(("trim", [
            "--sample-table", sample_table,
            "--adapters", cfg["adapters"],
            "--out-dir", str(out_dir / "trimmed"),
            "--trimmed-sample-table", str(trimmed_table),
            "--minimum-length", cfg.get("minimum_length", "50"),
            "--jobs", cfg.get("jobs", "4"),
            "--threads", cfg.get("threads", "2"),
        ], trimmed_table))
        align_table = trimmed_table
    else:
        align_table = Path(sample_table)

    plan.append(("align", [
        "--sample-table", str(align_table),
        "--peptide-table", peptide_table,
        "--out-dir", str(out_dir / "alignment"),
        "--mode", mode,
        "--preset", preset,
        "--jobs", cfg.get("jobs", "4"),
        "--threads", cfg.get("threads", "2"),
    ] + (["--scratch", cfg["scratch"]] if "scratch" in cfg else []), counts))

    plan.append(("score", [
        "--counts", str(counts),
        "--sample-table", str(align_table),
        "--peptide-table", peptide_table,
        "--out-dir", str(out_dir / "enrichment"),
        "--method", method,
        "--min-lib-size", cfg.get("min_lib_size", "500"),
        "--posterior-threshold", cfg.get("posterior_threshold", "0.5"),
        "--edger-fdr", cfg.get("edger_fdr", "0.05"),
        "--seed", cfg.get("seed", "20260101"),
    ] + (["--beads-rr"] if enabled(cfg, "beads_rr") else [])
      + (["--replicate-rule"] if enabled(cfg, "replicate_rule") else [])
      + [flag for key, name in (("threshold", "--threshold"),
                                ("gp_bins", "--gp-bins"),
                                ("gp_null", "--gp-null"),
                                ("replicate_group_column", "--replicate-group-column"))
         if key in cfg for flag in (name, cfg[key])], posterior))

    prioritise = [
        "--posterior", str(posterior),
        "--hits", str(hits),
        "--sample-table", str(align_table),
        "--peptide-table", peptide_table,
        "--out-dir", str(out_dir / "prioritised"),
        "--min-replicates", cfg.get("min_replicates", "2"),
        "--min-groups", cfg.get("min_groups", "2"),
        "--all-replicates", cfg.get("all_replicates", "3"),
        "--tile-step", cfg.get("tile_step", "28"),
    ]
    for key, flag in (("group_column", "--group-column"),
                      ("role_column", "--role-column"),
                      ("role", "--role"),
                      ("antigen_columns", "--antigen-columns"),
                      ("position_column", "--position-column")):
        if key in cfg:
            prioritise += [flag, cfg[key]]
    plan.append(("prioritise", prioritise,
                 out_dir / "prioritised" / "prioritisation_summary.json"))
    return plan


def adapter_lists(path):
    """Return the 5 prime adapter sequences per read, as the workflow expects.

    The workflow takes comma separated sequences rather than the stage route's
    CSV, and has no 3 prime parameter, so those rows are counted and dropped.
    """
    r1, r2, dropped = [], [], 0
    for row in csv.DictReader(open(path, newline="")):
        seq = (row.get("sequence") or "").strip().upper()
        if not seq:
            continue
        if (row.get("end") or "").strip().strip("'\"")[:1] != "5":
            dropped += 1
            continue
        read = (row.get("read") or "").strip().upper()
        if read == "R1":
            r1.append(seq)
        elif read == "R2":
            r2.append(seq)
    return r1, r2, dropped


def render_value(value):
    """Quote a value unless it is a Groovy boolean or number."""
    text = str(value)
    if text.lower() in ("true", "false"):
        return text.lower()
    try:
        float(text)
    except ValueError:
        return "'" + text.replace("\\", "\\\\").replace("'", "\\'") + "'"
    return text


def write_workflow_config(cfg, workflow, out_dir):
    """Render the [workflow] section as a Nextflow params file.

    Inputs come from [phipstream] so the two routes cannot disagree about which
    tables and adapters a dataset uses. Anything set in [workflow] wins.
    """
    params = {"sample_table": cfg["sample_table"],
              "peptide_table": cfg["peptide_table"],
              "results": str(out_dir / "workflow")}
    if "adapters" in cfg:
        r1, r2, dropped = adapter_lists(cfg["adapters"])
        if dropped:
            print(f"[stages] {dropped} three prime adapter row(s) dropped, "
                  "the workflow has no parameter for them")
        params.update(run_adapter_trim="true",
                      adapter_r1_5p_list=",".join(r1),
                      adapter_r2_5p_list=",".join(r2))
    if enabled(cfg, "fastqc"):
        params["fastqc_dir"] = str(out_dir / "fastqc")
    params.update(workflow)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "workflow.config"
    body = "\n".join(f"    {k} = {render_value(v)}"
                     for k, v in sorted(params.items()))
    path.write_text("params {\n" + body + "\n}\n")
    return path


def write_job_script(cfg, config_path, out_dir, project, python):
    name = f"phipstream_{Path(config_path).stem}_{datetime.now():%Y%m%d_%H%M%S}"
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    script = log_dir / f"{name}.qsub"
    modules = cfg.get("modules", "tools apptainer")
    body = f"""#!/usr/bin/env bash
#PBS -N {name}
#PBS -A {project}
#PBS -W group_list={project}
#PBS -l nodes=1:ppn={cfg.get('ppn', '8')},mem={cfg.get('mem', '32gb')},walltime={cfg.get('walltime', '24:00:00')}
#PBS -o {log_dir}/{name}.out
#PBS -e {log_dir}/{name}.err

set -euo pipefail
module load {modules}

cd {shlex.quote(str(BIN.parent))}
export PHIPSTREAM_PYTHON={shlex.quote(python)}
{shlex.quote(python)} {shlex.quote(str(BIN / 'run_stages.py'))} \\
    {shlex.quote(str(Path(config_path).resolve()))} --resume
"""
    script.write_text(body)
    script.chmod(0o755)
    return script


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="INI file with a [phipstream] section")
    ap.add_argument("--resume", action="store_true",
                    help="skip a stage whose output is already present")
    ap.add_argument("--submit", action="store_true",
                    help="write a PBS job script and submit it instead of running")
    ap.add_argument("--computerome_project", default="",
                    help="PBS account and group list, required with --submit")
    ap.add_argument("--workflow", action="store_true",
                    help="render the [workflow] section and run the bundled "
                         "workflow instead of the stages")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the stage commands without running them")
    args = ap.parse_args()

    cfg, workflow = load_config(args.config)
    out_dir = Path(require(cfg, "out_dir", args.config)).resolve()
    # A python key in the config wins, then the environment, then whatever is
    # running this script. The stage scripts need pandas, this one does not.
    python = cfg.get("python") or os.environ.get("PHIPSTREAM_PYTHON") or sys.executable

    if args.workflow:
        if not workflow:
            sys.exit(f"{args.config} has no [workflow] section")
        rendered = write_workflow_config(cfg, workflow, out_dir)
        print(f"[stages] wrote {rendered}")
        if args.dry_run:
            return 0
        argv = [str(BIN / "phipstream"),
                "submit" if args.submit else "nextflow", str(rendered)]
        if args.computerome_project:
            argv += ["--computerome_project", args.computerome_project]
        return subprocess.run(argv).returncode

    if args.submit:
        if not args.computerome_project:
            sys.exit("--submit needs --computerome_project")
        script = write_job_script(cfg, args.config, out_dir,
                                  args.computerome_project, python)
        print(f"[stages] wrote {script}")
        if subprocess.run(["which", "qsub"], capture_output=True).returncode == 0:
            job = subprocess.run(["qsub", str(script)], capture_output=True, text=True)
            if job.returncode != 0:
                sys.exit(f"qsub failed: {job.stderr.strip()}")
            print(f"[stages] submitted {job.stdout.strip()}")
        else:
            print("[stages] qsub not on PATH, submit the script from a node that has one")
        return 0

    plan = build_commands(cfg, args.config, out_dir)
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[stages] {len(plan)} stages into {out_dir}")
    for stage, options, marker in plan:
        argv = [python, str(BIN / SCRIPTS[stage])] + [str(x) for x in options]
        if args.dry_run:
            print(" ".join(shlex.quote(a) for a in argv))
            continue
        if args.resume and Path(marker).exists():
            print(f"[stages] {stage} already done, {marker} exists")
            continue
        print(f"[stages] {stage}")
        start = time.time()
        if subprocess.run(argv).returncode != 0:
            sys.exit(f"[stages] {stage} failed")
        print(f"[stages] {stage} finished in {time.time() - start:.0f}s")
    if not args.dry_run:
        print(f"[stages] complete, results under {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
