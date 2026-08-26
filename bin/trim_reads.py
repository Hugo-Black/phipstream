#!/usr/bin/env python3
"""Cutadapt wrapper for the trimming stage.

The command reads the shared sample table plus an adapter CSV, writes gzipped
FASTQ outputs, and emits a replacement sample table that points at those files.
The source table is read only.

Adapter CSV fields are name, read (R1/R2), end (5/3), sequence, and anchored.
  end 5 becomes cutadapt -g for R1 or -G for R2
  end 3 becomes cutadapt -a for R1 or -A for R2
  anchored true adds ^SEQ for 5 prime adapters or SEQ$ for 3 prime adapters
More than one row can fill the same adapter slot, which lets cutadapt choose the
best matching leader or tail for each read.
"""
import argparse
import csv
import shlex
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import containers  # noqa: E402

SLOTS = {("R1", "5"): "g", ("R1", "3"): "a", ("R2", "5"): "G", ("R2", "3"): "A"}


def load_adapters(path):
    slots = {"g": [], "a": [], "G": [], "A": []}
    for row in csv.DictReader(open(path, newline="")):
        seq = (row.get("sequence") or "").strip().upper()
        if not seq:
            continue
        read = (row.get("read") or "").strip().upper()
        end = (row.get("end") or "").strip().strip("'\"")[:1]
        name = (row.get("name") or "").strip() or seq[:10]
        anchored = (row.get("anchored") or "").strip().lower() in (
            "true", "t", "yes", "y", "1")
        key = SLOTS.get((read, end))
        if key is None:
            sys.exit(f"adapter row needs read R1/R2 and end 5/3, got {row!r}")
        if end == "5":
            spec = f"{name}=^{seq}" if anchored else f"{name}={seq}"
        else:
            spec = f"{name}={seq}$" if anchored else f"{name}={seq}"
        slots[key].append(spec)
    if not any(slots.values()):
        sys.exit(f"no adapters parsed from {path}")
    return slots


def build_command(slots, r1_in, r2_in, r1_out, r2_out, min_len, threads, json_path):
    parts = ["cutadapt"]
    for spec in slots["g"]:
        parts += ["-g", spec]
    for spec in slots["a"]:
        parts += ["-a", spec]
    if r2_in is not None:
        for spec in slots["G"]:
            parts += ["-G", spec]
        for spec in slots["A"]:
            parts += ["-A", spec]
    # cutadapt normally clips one adapter per pass. Add passes when both the
    # leader and the downstream adapter may be present on the same read.
    n_ends = bool(slots["g"] or slots["G"]) + bool(slots["a"] or slots["A"])
    if n_ends > 1:
        parts += ["-n", str(n_ends)]
    parts += ["--match-read-wildcards", "--minimum-length", str(min_len),
              "-j", str(threads), "--json", str(json_path), "-o", str(r1_out)]
    if r2_in is not None:
        parts += ["-p", str(r2_out), str(r1_in), str(r2_in)]
    else:
        parts += [str(r1_in)]
    return " ".join(shlex.quote(p) for p in parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-table", required=True,
                    help="table with fastq_filepath and optional fastq_r2_filepath")
    ap.add_argument("--adapters", required=True, help="adapter spec CSV")
    ap.add_argument("--out-dir", required=True, help="destination for trimmed reads")
    ap.add_argument("--trimmed-sample-table", required=True,
                    help="table written with fastq paths repointed at the trimmed reads")
    ap.add_argument("--minimum-length", type=int, default=50,
                    help="drop a read once trimming leaves it shorter than this. "
                         "It is the shortest insert still worth aligning, not a "
                         "fraction of the read length (default 50)")
    ap.add_argument("--jobs", type=int, default=4, help="samples trimmed at once")
    ap.add_argument("--threads", type=int, default=2, help="cutadapt threads per sample")
    args = ap.parse_args()

    slots = load_adapters(args.adapters)
    out_dir = Path(args.out_dir).resolve()
    log_dir = out_dir / "cutadapt_logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.sample_table)))
    if not rows:
        sys.exit("empty sample table")
    fields = list(rows[0].keys())
    if "fastq_filepath" not in fields:
        sys.exit("sample table needs a fastq_filepath column")

    jobs, binds = [], {str(out_dir)}
    for row in rows:
        r1_in = Path(row["fastq_filepath"]).resolve()
        r2_raw = (row.get("fastq_r2_filepath") or "").strip()
        r2_in = Path(r2_raw).resolve() if r2_raw else None
        for path in (r1_in, r2_in):
            if path is not None:
                if not path.exists():
                    sys.exit(f"missing read file: {path}")
                binds.add(str(path.parent))
        r1_out = out_dir / r1_in.name
        r2_out = out_dir / r2_in.name if r2_in is not None else None
        stem = r1_in.name.replace(".fastq.gz", "").replace(".fq.gz", "")
        cmd = build_command(slots, r1_in, r2_in, r1_out, r2_out,
                            args.minimum_length, args.threads,
                            log_dir / f"{stem}.cutadapt.json")
        jobs.append((row, cmd, log_dir / f"{stem}.cutadapt.log", r1_out, r2_out))

    n_pe = sum(1 for j in jobs if j[4] is not None)
    print(f"[trim] {len(jobs)} samples ({n_pe} paired, {len(jobs) - n_pe} single) "
          f"into {out_dir}")
    print(f"[trim] adapters: -g {len(slots['g'])} -a {len(slots['a'])} "
          f"-G {len(slots['G'])} -A {len(slots['A'])}")

    def trim_one(job):
        _row, cmd, log, _r1, _r2 = job
        proc = containers.run("cutadapt", cmd, binds, out_dir,
                              capture=True, check=False)
        log.write_text((proc.stdout or "") + (proc.stderr or ""))
        return proc.returncode

    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(trim_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futures), 1):
            job = futures[fut]
            name = Path(job[0]["fastq_filepath"]).name
            code = fut.result()
            if code != 0:
                failures.append(name)
            print(f"  [{i}/{len(jobs)}] {name}: "
                  f"{'ok' if code == 0 else f'failed with exit code {code}'}")
    if failures:
        sys.exit(f"{len(failures)} sample(s) failed, see {log_dir}")

    for row, _cmd, _log, r1_out, r2_out in jobs:
        row["fastq_filepath"] = str(r1_out)
        if r2_out is not None and "fastq_r2_filepath" in fields:
            row["fastq_r2_filepath"] = str(r2_out)
    out_table = Path(args.trimmed_sample_table)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    with open(out_table, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[trim] wrote {out_table}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
