#!/usr/bin/env python3
"""Create lightweight hit tables for the quick end to end test path.

The rule labels peptides whose counts sit far above the sample median, which is
enough to recover the deliberately boosted fixture tiles. This is only a test
shortcut for exercising prioritisation. Set BEER=1 to run the real scorer.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--sample-table", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fold", type=float, default=8.0,
                    help="multiple of the sample median that counts as enriched")
    args = ap.parse_args()

    counts = pd.read_csv(args.counts, index_col=0)
    table = pd.read_csv(args.sample_table, dtype=str)
    key = dict(zip(table["sample_name"], table["technical_replicate_id"]))
    beads = {key[n] for n, s in zip(table["sample_name"], table.get(
        "control_status", pd.Series("", index=table.index)))
        if s == "beads_only" and n in key}

    counts = counts.rename(columns=key)
    threshold = counts.median(axis=0) * args.fold
    hits = counts.gt(threshold, axis=1)
    posterior = hits.astype(float).mask(hits, 0.99).mask(~hits, 0.01)
    for column in counts.columns:
        if column in beads:
            posterior[column] = float("nan")
            hits[column] = False

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    posterior.index.name = "peptide_id"
    hits.index.name = "peptide_id"
    posterior.to_csv(out / "beer_posterior.csv.gz")
    hits.astype(int).to_csv(out / "beer_hits.csv.gz")
    print(f"[test] stand-in calls: {int(hits.values.sum())} over "
          f"{len(counts.columns) - len(beads)} replicates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
