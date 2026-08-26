#!/usr/bin/env python3
"""Write an nf-core/seqinspector input sheet from the project sample table.

Seqinspector is launched outside this script and expects its own CSV layout.
This converter keeps QC aligned with the pipeline by drawing read paths and
labels from the same sample metadata used by later stages.

The output columns are sample, fastq_1, fastq_2, rundir, and tags. Tags carry a
sample class or control label for grouped reporting.
"""
import argparse
import csv
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-table", required=True,
                    help="table with fastq_filepath and a name column")
    ap.add_argument("--output", required=True,
                    help="sample sheet to write")
    ap.add_argument("--reads-prefix", default="",
                    help="prefix joined to relative fastq paths")
    ap.add_argument("--tags-column", default="control_status",
                    help="column used as the grouping tag (default control_status)")
    ap.add_argument("--name-column", default="sample_name",
                    help="column used as the sample name (default sample_name)")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.sample_table)))
    if not rows:
        sys.exit("empty sample table")
    if "fastq_filepath" not in rows[0]:
        sys.exit("sample table needs a fastq_filepath column")
    name_column = (args.name_column if args.name_column in rows[0]
                   else "technical_replicate_id")

    def resolve(value):
        value = (value or "").strip()
        if not value:
            return ""
        path = Path(value)
        if args.reads_prefix and not path.is_absolute():
            path = Path(args.reads_prefix) / path
        return str(path.resolve())

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", "fastq_1", "fastq_2", "rundir", "tags"])
        writer.writeheader()
        for i, row in enumerate(rows, 1):
            fastq_1 = resolve(row.get("fastq_filepath"))
            if not fastq_1:
                continue
            writer.writerow({
                "sample": (row.get(name_column) or f"sample_{i}").strip(),
                "fastq_1": fastq_1,
                "fastq_2": resolve(row.get("fastq_r2_filepath")),
                "rundir": "",
                "tags": (row.get(args.tags_column) or "").strip(),
            })
            written += 1
    print(f"[qc] wrote {out} ({written} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
