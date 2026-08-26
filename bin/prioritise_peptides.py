#!/usr/bin/env python3
"""Build replicate concordance summaries and peptide shortlists.

The first pass measures concordance inside each replicate group that has enough
fitted replicates. It asks what fraction of peptides called anywhere in the
group are supported by the requested number of replicates. Groups without calls
are omitted from the concordance table.

The second pass starts from peptides with enough replicate support in at least
one group. It keeps peptides that recur across groups, appear in every fitted
replicate of a group, or have a half-overlapping tile called in the same group.
External labels are not used by these rules.

Replicates without posterior values are ignored, because the hit table encodes
missing model output and true negative calls in the same way.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def load_calls(posterior_path, hits_path, threshold):
    """Return the boolean hit matrix and columns with fitted posterior values.

    BEER can omit posterior values for super-enriched peptides. Those cells are
    treated as certainty before fitted replicate columns are detected, preventing
    high-confidence calls from being lost.
    """
    posterior = pd.read_csv(posterior_path, index_col=0)
    hits = pd.read_csv(hits_path, index_col=0)
    hits = hits.reindex(index=posterior.index, columns=posterior.columns)
    hits = hits.fillna(0).astype(float) > 0
    filled = posterior.mask(posterior.isna() & hits, 1.0)
    modelled = [c for c in filled.columns if filled[c].notna().any()]
    if threshold is not None:
        hits = hits | (filled > threshold).fillna(False)
    return hits, modelled


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--posterior", required=True, help="beer_posterior.csv.gz")
    ap.add_argument("--hits", required=True, help="beer_hits.csv.gz")
    ap.add_argument("--sample-table", required=True,
                    help="sample table with technical_replicate_id")
    ap.add_argument("--peptide-table", required=True,
                    help="peptide table with peptide_id and the grouping columns")
    ap.add_argument("--out-dir", required=True,
                    help="destination for the concordance and shortlist tables")
    ap.add_argument("--group-column", default="participant_ID",
                    help="column that names the replicate group (default participant_ID)")
    ap.add_argument("--role-column", default="sample_role",
                    help="column that names the sample class (default sample_role)")
    ap.add_argument("--role", default="serum",
                    help="sample class to prioritise, empty for all (default serum)")
    ap.add_argument("--min-replicates", type=int, default=2,
                    help="replicates of one group a peptide must be enriched in (default 2)")
    ap.add_argument("--min-groups", type=int, default=2,
                    help="groups a peptide must recur across to be kept (default 2)")
    ap.add_argument("--all-replicates", type=int, default=3,
                    help="replicates that count as enrichment in every replicate (default 3)")
    ap.add_argument("--tile-step", type=int, default=28,
                    help="tile start offset of an adjacent overlapping tile (default 28)")
    ap.add_argument("--antigen-columns", default="virus,antigen",
                    help="peptide table columns an adjacent tile must share")
    ap.add_argument("--position-column", default="tile_start",
                    help="peptide table column holding the tile start")
    ap.add_argument("--posterior-threshold", type=float, default=None,
                    help="recompute calls at this posterior instead of trusting the hit matrix")
    args = ap.parse_args()

    hits, modelled = load_calls(args.posterior, args.hits, args.posterior_threshold)
    table = pd.read_csv(args.sample_table, dtype=str)
    if "technical_replicate_id" not in table.columns:
        sys.exit("sample table needs a technical_replicate_id column")
    table["replicate_id"] = table["technical_replicate_id"].str.strip()
    hits.columns = [str(c).strip() for c in hits.columns]
    modelled = {str(c).strip() for c in modelled}

    peptides = pd.read_csv(args.peptide_table).set_index("peptide_id")
    keys = [c.strip() for c in args.antigen_columns.split(",") if c.strip()]
    for column in keys + [args.position_column]:
        if column not in peptides.columns:
            sys.exit(f"peptide table has no {column} column")

    selected = table
    if args.role:
        if args.role_column not in table.columns:
            sys.exit(f"--role {args.role!r} needs a {args.role_column} column in "
                     "the sample table. Pass --role '' to use every sample")
        selected = table[table[args.role_column].fillna("").str.strip() == args.role]
    if args.group_column not in selected.columns:
        sys.exit(f"sample table has no {args.group_column} column")
    selected = selected[selected[args.group_column].fillna("").str.strip() != ""]
    if selected.empty:
        sys.exit(f"no samples matched {args.role_column}={args.role!r}")

    groups = {}
    for name, rows in selected.groupby(args.group_column):
        reps = [r for r in rows["replicate_id"]
                if r in modelled and r in hits.columns]
        groups[name] = reps

    all_reps = [r for reps in groups.values() for r in reps]
    n_events = int(hits[all_reps].sum().sum())
    detected = sorted(hits.index[hits[all_reps].any(axis=1)])

    concordance, supported = [], []
    for name in sorted(groups):
        reps = groups[name]
        row = {"group": name, "n_samples": int((selected[args.group_column] == name).sum()),
               "n_modelled": len(reps), "n_union": "", "n_repeated": "", "concordance": ""}
        if len(reps) >= args.min_replicates:
            support = hits[reps].sum(axis=1)
            union = int((support >= 1).sum())
            repeated = int((support >= args.min_replicates).sum())
            row.update(n_union=union, n_repeated=repeated,
                       concordance=(repeated / union) if union else "")
            for pid in support.index[support >= args.min_replicates]:
                supported.append({"peptide_id": pid, "group": name,
                                  "support": int(support[pid])})
        concordance.append(row)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    conc = pd.DataFrame(concordance)
    conc.to_csv(out_dir / "replicate_concordance.csv", index=False)
    cons = pd.DataFrame(supported, columns=["peptide_id", "group", "support"])
    cons.to_csv(out_dir / "replicate_supported.csv", index=False)

    if cons.empty:
        print("[prioritise] no peptide reached the replicate threshold")
        criteria, pairs = {"groups": set(), "replicates": set(), "adjacent": set()}, []
    else:
        cons = cons.join(peptides[keys + [args.position_column]], on="peptide_id")
        by_group = cons.groupby("peptide_id")["group"].nunique()
        c_groups = set(by_group.index[by_group >= args.min_groups])
        c_reps = set(cons.loc[cons["support"] >= args.all_replicates, "peptide_id"])
        c_adj, pair_set = set(), set()
        for _name, block in cons.groupby("group"):
            for row in block.itertuples():
                same = block
                for column in keys:
                    same = same[same[column] == getattr(row, column)]
                offset = (same[args.position_column]
                          - getattr(row, args.position_column)).abs()
                for other in same.loc[offset == args.tile_step, "peptide_id"]:
                    c_adj.update((row.peptide_id, other))
                    pair_set.add(tuple(sorted((row.peptide_id, other))))
        criteria = {"groups": c_groups, "replicates": c_reps, "adjacent": c_adj}
        pairs = sorted(pair_set)

    shortlist = sorted(set().union(*criteria.values())) if cons.shape[0] else []
    columns = [c for c in ("virus", "antigen", "region", args.position_column)
               if c in peptides.columns]
    rows = []
    for pid in shortlist:
        block = cons[cons["peptide_id"] == pid]
        rows.append({
            "peptide_id": pid,
            **{c: peptides.at[pid, c] for c in columns},
            "groups": ",".join(sorted(block["group"])),
            "max_support": int(block["support"].max()),
            "criterion_groups": pid in criteria["groups"],
            "criterion_all_replicates": pid in criteria["replicates"],
            "criterion_adjacent_tile": pid in criteria["adjacent"],
        })
    pd.DataFrame(rows).to_csv(out_dir / "shortlist.csv", index=False)
    pd.DataFrame(pairs, columns=["peptide_id_a", "peptide_id_b"]).to_csv(
        out_dir / "adjacent_pairs.csv", index=False)

    finite = pd.to_numeric(conc["concordance"], errors="coerce").dropna()
    summary = {
        "role": args.role,
        "n_groups": len(groups),
        "n_groups_scored": int(len(finite)),
        "n_modelled_replicates": len(all_reps),
        "n_calls": n_events,
        "n_peptides_detected": len(detected),
        "n_peptides_replicate_supported": int(cons["peptide_id"].nunique()) if len(cons) else 0,
        "n_shortlist": len(shortlist),
        "criteria": {k: len(v) for k, v in criteria.items()},
        "concordance_median": float(finite.median()) if len(finite) else None,
        "settings": {"min_replicates": args.min_replicates,
                     "min_groups": args.min_groups,
                     "all_replicates": args.all_replicates,
                     "tile_step": args.tile_step},
    }
    (out_dir / "prioritisation_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[prioritise] {summary['n_calls']} calls over "
          f"{summary['n_modelled_replicates']} modelled replicates")
    print(f"[prioritise] {summary['n_peptides_detected']} peptides detected, "
          f"{summary['n_peptides_replicate_supported']} replicate supported, "
          f"{summary['n_shortlist']} shortlisted")
    print(f"[prioritise] concordance median "
          f"{summary['concordance_median']} over {summary['n_groups_scored']} groups")
    return 0


if __name__ == "__main__":
    sys.exit(main())
