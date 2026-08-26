#!/usr/bin/env python3
"""FastQC wrapper for the read quality stage.

The command reads the shared sample table and writes one FastQC archive per read
file. The contamination diagnostic in the bundled workflow reads those archives
rather than the reads themselves, so this stage has to run before it. Point the
workflow's fastqc_dir at this stage's output directory.

Run it on untrimmed reads. Trimming removes the leader that overrepresented
sequence reporting is most likely to surface.
"""
import argparse
import csv
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import containers  # noqa: E402

SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def read_paths(table):
    """Return the read files named by the sample table, R1 before R2."""
    rows = list(csv.DictReader(open(table)))
    if not rows:
        sys.exit("empty sample table")
    if "fastq_filepath" not in rows[0]:
        sys.exit("sample table needs a fastq_filepath column")
    paths = []
    for row in rows:
        for key in ("fastq_filepath", "fastq_r2_filepath"):
            raw = (row.get(key) or "").strip()
            if not raw:
                continue
            path = Path(raw).resolve()
            if not path.exists():
                sys.exit(f"missing read file: {path}")
            if path not in paths:
                paths.append(path)
    return paths


def archive_for(path, out_dir):
    name = path.name
    for suffix in SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return out_dir / f"{name}_fastqc.zip"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-table", required=True,
                    help="table with fastq_filepath and optional fastq_r2_filepath")
    ap.add_argument("--out-dir", required=True,
                    help="destination for the FastQC archives")
    ap.add_argument("--threads", type=int, default=4,
                    help="read files processed at once (default 4)")
    args = ap.parse_args()

    paths = read_paths(args.sample_table)
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    binds = {str(out_dir)} | {str(p.parent) for p in paths}

    print(f"[fastqc] {len(paths)} read files into {out_dir}")
    argv = ["fastqc", "-t", str(args.threads), "-o", str(out_dir)]
    argv += [str(p) for p in paths]
    containers.run("fastqc", " ".join(shlex.quote(a) for a in argv), binds, out_dir)

    index = out_dir / "fastqc_files.csv"
    missing = []
    with open(index, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["read_file", "archive"])
        for path in paths:
            archive = archive_for(path, out_dir)
            if not archive.exists():
                missing.append(path.name)
            writer.writerow([str(path), str(archive)])
    if missing:
        sys.exit(f"{len(missing)} archive(s) not produced, first is {missing[0]}")
    print(f"[fastqc] wrote {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
