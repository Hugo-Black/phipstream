#!/usr/bin/env python3
"""Run the enrichment scoring stage through edgeR and BEER.

The wrapper prepares the three wide tables expected by the R scripts, fits the
beads-only background with edgeR, serializes the PhIPData object, and then sends
that object through BEER. Replicate calls come from the posterior threshold, with
BEER super-enriched records retained when a posterior value is absent.

Samples below --min-lib-size are excluded before fitting and restored later as
empty columns. That preserves the original sample shape while keeping shallow
samples out of the model. Unfitted replicates stay as missing posterior columns
rather than being converted to zeroes.

A fixed --seed is supplied unless the caller overrides it, and the selected value
is written into the run summary.
"""
import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import containers  # noqa: E402

R_DIR = Path(__file__).resolve().parent.parent / "phip-flow" / "bin"
# Score matrix and hit matrix each method writes, and what prioritise reads.
SCORE_FILES = {"beer": ("beer_posterior.csv.gz", "beer_hits.csv.gz"),
               "edger": ("edger_logpval.csv.gz", "edger_hits.csv.gz")}
TRUE_WORDS = ("true", "1", "yes", "t")


def is_true(value):
    return str(value).strip().lower() in TRUE_WORDS


def ensure_control_status(table):
    """Return a table with a usable control_status column.

    Existing non-empty status values are preserved. Otherwise the status is
    inferred from role and boolean metadata columns.
    """
    df = table.copy()
    if ("control_status" in df.columns
            and df["control_status"].fillna("").str.strip().ne("").any()):
        return df
    idx = df.index
    role = (df["sample_role"].fillna("").str.lower()
            if "sample_role" in df.columns else pd.Series("", index=idx))
    background = (df["in_background_null"].map(is_true)
                  if "in_background_null" in df.columns
                  else pd.Series(False, index=idx))
    input_ref = (df["is_input_reference"].map(is_true)
                 if "is_input_reference" in df.columns
                 else pd.Series(False, index=idx))
    status = pd.Series("empirical", index=idx)
    status[role.str.contains("library") | input_ref] = "library"
    status[background | role.str.contains("beads")] = "beads_only"
    status[role.eq("undetermined")] = "undetermined"
    df["control_status"] = status.values
    print(f"[score] derived control_status: {status.value_counts().to_dict()}")
    return df


def load_inputs(counts_path, sample_table_path, peptide_table_path):
    """Load counts and metadata, then align them on integer sample identifiers."""
    counts = pd.read_csv(counts_path)
    counts = counts.rename(columns={counts.columns[0]: "peptide_id"})
    counts["peptide_id"] = counts["peptide_id"].astype(int)
    counts = counts.set_index("peptide_id")

    table = ensure_control_status(pd.read_csv(sample_table_path, dtype=str))
    for column in ("sample_role", "control_status"):
        if column in table.columns:
            table = table[table[column].fillna("").str.strip().str.lower()
                          != "undetermined"]
    table = table.reset_index(drop=True)
    if "sample_name" not in table.columns:
        sys.exit("sample table needs a sample_name column")

    key = "technical_replicate_id"
    usable = (key in table.columns
              and table[key].notna().all()
              and table[key].str.strip().str.lstrip("-").str.isdigit().all()
              and table[key].astype(int).is_unique)
    table["sample_id"] = (table[key].astype(int) if usable
                          else range(1, len(table) + 1))
    print(f"[score] sample_id from "
          f"{'technical_replicate_id' if usable else 'row order'}")

    name_to_id = dict(zip(table["sample_name"], table["sample_id"]))
    keep = [c for c in counts.columns if c in name_to_id]
    dropped = [c for c in counts.columns if c not in name_to_id]
    if dropped:
        print(f"[score] {len(dropped)} count column(s) absent from the sample "
              f"table, dropped: {dropped[:8]}")
    counts = counts[keep].rename(columns=name_to_id)
    counts.columns = counts.columns.astype(int)

    peptides = pd.read_csv(peptide_table_path)
    peptides["peptide_id"] = peptides["peptide_id"].astype(int)
    peptides = peptides.set_index("peptide_id")

    counts = counts.reindex(peptides.index).fillna(0).astype(int).sort_index()
    counts = counts[sorted(counts.columns)]
    samples = table[table["sample_name"].isin(keep)].set_index("sample_id")
    samples = samples.sort_index()
    return counts, samples, peptides.loc[counts.index]


def write_r_inputs(work, counts, samples, peptides):
    """Write the count, sample, and peptide tables consumed by the R scripts."""
    counts_out = counts.copy()
    counts_out.index.name = ""
    counts_out.to_csv(work / "dataset_counts.csv", na_rep="NA")
    samples.convert_dtypes().to_csv(work / "dataset_sample_annotation_table.csv",
                                    index_label="sample_id")
    peptides.convert_dtypes().to_csv(work / "dataset_peptide_annotation_table.csv",
                                     index_label="peptide_id")


def restore_shape(frame, columns):
    """Put an R result back onto the complete set of integer sample IDs."""
    frame = frame.copy()
    frame.index = frame.index.astype(int)
    frame.columns = [int(str(c).lstrip("X")) for c in frame.columns]
    return frame.reindex(columns=columns)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="peptide by sample counts CSV")
    ap.add_argument("--sample-table", required=True,
                    help="sample table carrying control_status or the columns "
                         "it is derived from")
    ap.add_argument("--peptide-table", required=True,
                    help="peptide table with a peptide_id column")
    ap.add_argument("--out-dir", required=True,
                    help="destination for the score and hit matrices")
    ap.add_argument("--method", choices=("beer", "edger"), default="beer",
                    help="scoring method. beer runs edgeR to fit the prior and "
                         "then BEER over it. edger stops after edgeR, which is "
                         "far quicker and suits a deployment check (default beer)")
    ap.add_argument("--min-lib-size", type=int, default=500,
                    help="drop a sample below this many library-mapped reads (default 500)")
    ap.add_argument("--posterior-threshold", type=float, default=0.5,
                    help="posterior probability above which a peptide is called (default 0.5)")
    ap.add_argument("--beads-rr", action="store_true",
                    help="run each beads-only sample against the others")
    ap.add_argument("--edger-fdr", type=float, default=0.05,
                    help="BH cutoff for the edgeR hit matrix (default 0.05)")
    ap.add_argument("--seed", type=int, default=20260101,
                    help="sampler seed, pinned so a rerun reproduces the calls. "
                         "The default value carries no meaning beyond being fixed")
    ap.add_argument("--keep-workdir", default="",
                    help="copy the R working directory here for inspection")
    args = ap.parse_args()

    counts, samples, peptides = load_inputs(args.counts, args.sample_table,
                                            args.peptide_table)
    n_beads = int((samples["control_status"] == "beads_only").sum())
    print(f"[score] counts {counts.shape[0]} peptides by {counts.shape[1]} samples, "
          f"{n_beads} beads-only")
    if n_beads < 2:
        sys.exit(f"BEER needs a beads-only group, found {n_beads} sample(s). "
                 "Mark them with control_status beads_only")
    if n_beads < 4:
        print(f"[score] warning: {n_beads} beads-only samples. Four to eight mock "
              "immunoprecipitations give a more reliable prior", file=sys.stderr)

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="phipstream_beer_"))
    try:
        write_r_inputs(work, counts, samples, peptides)
        needed = ("run_edgeR.Rscript", "run_BEER.Rscript") if args.method == "beer" \
            else ("run_edgeR.Rscript",)
        for script in needed:
            source = R_DIR / script
            if not source.exists():
                sys.exit(f"missing R step: {source}")
            shutil.copy(source, work / script)
        command = f"Rscript run_edgeR.Rscript {args.edger_fdr}"
        if args.method == "beer":
            command += (f" && Rscript run_BEER.Rscript {args.min_lib_size} "
                        f"{args.posterior_threshold} "
                        f"{'TRUE' if args.beads_rr else 'FALSE'} {args.seed}")
            print(f"[score] running edgeR then BEER, seed {args.seed}")
        else:
            print(f"[score] running edgeR only, BH cutoff {args.edger_fdr}")
        log = containers.run("beer", command, {str(work)}, work,
                             capture=True, check=False)
        (out_dir / "scoring.log").write_text((log.stdout or "") + (log.stderr or ""))
        if log.returncode != 0:
            sys.stderr.write((log.stdout or "")[-3000:])
            sys.stderr.write((log.stderr or "")[-3000:])
            sys.exit(f"{args.method} step failed, full log in "
                     f"{out_dir / 'scoring.log'}")

        columns = list(counts.columns)
        written = {}
        for name, target in (("beer_prob", "beer_posterior"),
                             ("beer_hits", "beer_hits"),
                             ("beer_fc_cond", "beer_fold_change_conditional"),
                             ("beer_fc_marg", "beer_fold_change_marginal"),
                             ("edgeR_hits", "edger_hits"),
                             ("edgeR_logfc", "edger_logfc"),
                             ("edgeR_logpval", "edger_logpval")):
            path = work / f"{name}.csv"
            if not path.exists():
                continue
            frame = restore_shape(pd.read_csv(path, index_col=0), columns)
            frame.index.name = "peptide_id"
            frame.to_csv(out_dir / f"{target}.csv.gz")
            written[target] = list(frame.shape)
        if args.keep_workdir:
            shutil.copytree(work, args.keep_workdir, dirs_exist_ok=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    score_name, hits_name = SCORE_FILES[args.method]
    posterior = pd.read_csv(out_dir / score_name, index_col=0)
    hits = pd.read_csv(out_dir / hits_name, index_col=0)
    modelled = [c for c in posterior.columns if posterior[c].notna().any()]
    depth = counts.sum(axis=0)
    summary = {
        "counts_source": str(Path(args.counts).resolve()),
        "sample_table": str(Path(args.sample_table).resolve()),
        "n_peptides": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "n_beads_only": n_beads,
        "method": args.method,
        "score_matrix": score_name,
        "hits_matrix": hits_name,
        "min_lib_size": args.min_lib_size,
        "posterior_threshold": args.posterior_threshold,
        "beads_rr": bool(args.beads_rr),
        "seed": args.seed,
        "images": containers.IMAGES,
        "n_below_min_lib_size": int((depth < args.min_lib_size).sum()),
        "n_modelled": len(modelled),
        "n_unmodelled": int(posterior.shape[1] - len(modelled)),
        "n_calls": int((hits.fillna(0) > 0).sum().sum()),
        "outputs": written,
    }
    (out_dir / "enrichment_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[score] {summary['n_modelled']} of {summary['n_samples']} replicates "
          f"modelled, {summary['n_calls']} calls")
    print(f"[score] wrote {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
