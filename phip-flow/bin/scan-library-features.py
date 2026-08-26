#!/usr/bin/env python3
"""Calculate sequence QC features for a peptide library.

The output table marks oligos whose composition suggests possible synthesis
problems, allowing review before ordering the library. A predicted-fail FASTA
and summary JSON are written beside the feature table.

Usage:
    scan-library-features.py \n        --peptide-table <path>/peptide_table.csv \n        --peptide-seq-col oligo \n        --peptide-protein-col peptide_translate \n        --peptide-id-col peptide_id \n        --trim5 26 --trim3 28 \n        --out-dir /tmp/lib_qc

When --peptide-table is a FASTA file, headers are used as ids.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

# Make the sibling helper module importable.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from library_qc_lib import (  # noqa: E402
    DEFAULT_THRESHOLDS,
    featurize_oligo,
    parse_fasta,
    predicted_fail,
)


def _read_peptide_table(
    path: Path,
    seq_col: str,
    protein_col: str | None,
    id_col: str | None,
    organism_col: str | None,
) -> list[dict]:
    """Load CSV or FASTA input as a list of peptide records."""
    if path.suffix.lower() in {".fasta", ".fa", ".fna"}:
        fasta = parse_fasta(path)
        return [
            {
                "peptide_id": header,
                seq_col: seq,
            }
            for header, seq in fasta.items()
        ]

    rows: list[dict] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if seq_col not in (reader.fieldnames or []):
            avail = ", ".join(reader.fieldnames or [])
            raise SystemExit(
                f"ERROR: --peptide-seq-col '{seq_col}' not in peptide table. "
                f"Available columns: {avail}"
            )
        for r in reader:
            rows.append(r)
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peptide-table", type=Path, required=True)
    p.add_argument("--peptide-seq-col", default="oligo",
                   help="DNA column name (default: oligo)")
    p.add_argument("--peptide-protein-col", default="peptide_translate",
                   help="Protein column name (default: peptide_translate); "
                        "if absent the script translates from the DNA column")
    p.add_argument("--peptide-id-col", default="peptide_id",
                   help="Peptide identifier column (default: peptide_id)")
    p.add_argument("--peptide-organism-col", default="virus",
                   help="Optional organism / virus column for per-group report panels")
    p.add_argument("--trim5", type=int, default=0,
                   help="Strip this many bp from the 5' end before scoring")
    p.add_argument("--trim3", type=int, default=0,
                   help="Strip this many bp from the 3' end before scoring")
    p.add_argument("--gc-window", type=int, default=30,
                   help="Sliding-window width for GC scan (default: 30)")
    p.add_argument("--thresholds", default="",
                   help="JSON map of threshold overrides")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    thresholds = dict(DEFAULT_THRESHOLDS)
    if args.thresholds:
        thresholds.update(json.loads(args.thresholds))

    rows = _read_peptide_table(
        args.peptide_table,
        args.peptide_seq_col,
        args.peptide_protein_col,
        args.peptide_id_col,
        args.peptide_organism_col,
    )

    if not rows:
        raise SystemExit(f"ERROR: peptide table at {args.peptide_table} is empty")

    feature_rows: list[dict] = []
    n_flagged = 0
    reason_counts: dict[str, int] = {}

    for r in rows:
        dna = (r.get(args.peptide_seq_col) or "").strip()
        if not dna:
            continue
        protein = (r.get(args.peptide_protein_col) or "").strip() or None
        feats = featurize_oligo(
            dna,
            protein=protein,
            trim5=args.trim5,
            trim3=args.trim3,
            gc_window=args.gc_window,
        )
        pf, reasons = predicted_fail(feats, thresholds)
        row_out = {
            "peptide_id": str(r.get(args.peptide_id_col, "")) or str(len(feature_rows)),
            "organism": str(r.get(args.peptide_organism_col, "")),
            **{k: v for k, v in feats.items() if k != "protein"},
            "predicted_fail": int(pf),
            "predicted_fail_reasons": ";".join(reasons),
            "protein": feats["protein"],
            "oligo_full": dna,
        }
        feature_rows.append(row_out)
        if pf:
            n_flagged += 1
            for reason in reasons:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1

    # Write the feature table.
    feat_path = args.out_dir / "library_features.csv"
    columns = [
        "peptide_id", "organism", "dna_len", "protein_len",
        "gc_overall", "gc_window_max", "gc_window_min", "gc_window_mean",
        "homopolymer_max_A", "homopolymer_max_C", "homopolymer_max_G",
        "homopolymer_max_T", "homopolymer_max",
        "dinucleotide_entropy_min", "tandem_repeat_span_dna",
        "aa_repeat_span", "aa_repeat_class",
        "predicted_fail", "predicted_fail_reasons",
        "protein", "oligo_full",
    ]
    with open(feat_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        for row in feature_rows:
            w.writerow({k: row.get(k, "") for k in columns})

    # Write predicted-fail oligos as FASTA records.
    fail_fa = args.out_dir / "library_predicted_fail.fasta"
    with open(fail_fa, "w") as fh:
        for row in feature_rows:
            if row["predicted_fail"]:
                reasons = row["predicted_fail_reasons"]
                fh.write(f">{row['peptide_id']} reasons={reasons} aa_repeat={row['aa_repeat_class']}\n")
                fh.write(f"{row['oligo_full']}\n")

    # Write top-level summary metrics.
    n_total = len(feature_rows)
    summary = {
        "n_peptides_scanned": n_total,
        "n_predicted_fail": n_flagged,
        "pct_predicted_fail": (100.0 * n_flagged / n_total) if n_total else 0.0,
        "fail_reasons_count": reason_counts,
        "thresholds_used": thresholds,
        "input_peptide_table": str(args.peptide_table.resolve()),
        "trim5": args.trim5,
        "trim3": args.trim3,
    }
    if n_total:
        # Add percentile summaries for numeric feature columns.
        def _pct(values: list[float], q: float) -> float:
            if not values:
                return 0.0
            vs = sorted(values)
            k = int(round(q * (len(vs) - 1)))
            return vs[k]

        for col in (
            "gc_overall", "gc_window_max", "homopolymer_max",
            "dinucleotide_entropy_min", "tandem_repeat_span_dna",
            "aa_repeat_span",
        ):
            vals = [float(row[col]) for row in feature_rows if row[col] not in ("", None)]
            summary[f"{col}_p05"] = _pct(vals, 0.05)
            summary[f"{col}_p50"] = _pct(vals, 0.50)
            summary[f"{col}_p95"] = _pct(vals, 0.95)

    sum_path = args.out_dir / "library_features_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2))

    # Print a concise terminal summary.
    print(f"Scanned {n_total} peptides into {feat_path}")
    print(f"Predicted-fail: {n_flagged} ({summary['pct_predicted_fail']:.2f}%)")
    if reason_counts:
        print("Per-reason counts:")
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"  {reason}: {count}")
    print(f"Predicted-fail FASTA: {fail_fa}")
    print(f"Summary JSON: {sum_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
