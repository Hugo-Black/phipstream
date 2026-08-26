#!/usr/bin/env python3
"""Generate the synthetic inputs used by the end to end test.

A fixed seed creates all reads and library rows on demand, which keeps binary
fixtures out of the repository. The layout matches the pipeline contract rather
than a biological study: antigens are tiled at a fixed offset, several groups
share one boosted peptide, and the first group includes an overlapping pair.
"""
import argparse
import csv
import gzip
import random
from pathlib import Path

ADAPTER = "TTGTTATTACTCGCGGCCCACTGCAG"
OLIGO_LEN = 168
TILE_STEP = 28
READ_LEN = 151
BASES = "ACGT"


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def build_library(rng, n_oligos, n_antigens):
    """Create library rows whose tile starts advance by a fixed offset."""
    rows = []
    per = max(1, n_oligos // n_antigens)
    for i in range(n_oligos):
        antigen = f"AG{i // per + 1:02d}"
        index = i % per
        rows.append({
            "peptide_id": i,
            "oligo": "".join(rng.choice(BASES) for _ in range(OLIGO_LEN)),
            "virus": f"V{(i // per) % 2 + 1}",
            "antigen": antigen,
            "region": f"{index * TILE_STEP}-{index * TILE_STEP + 56}",
            "tile_start": index * TILE_STEP,
        })
    return rows


def write_reads(path, oligos, counts, rng, error_rate=0.002, n_rate=0.004):
    """Write synthetic reads containing the leader followed by oligo sequence."""
    with gzip.open(path, "wt") as handle:
        read = 0
        for pid, n in counts.items():
            template = ADAPTER + oligos[pid]
            for _ in range(n):
                seq = list(template[:READ_LEN])
                for i, base in enumerate(seq):
                    draw = rng.random()
                    if draw < n_rate:
                        seq[i] = "N"
                    elif draw < n_rate + error_rate:
                        seq[i] = rng.choice(BASES.replace(base, ""))
                read += 1
                handle.write(f"@read{read}\n{''.join(seq)}\n+\n{'I' * len(seq)}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--oligos", type=int, default=200)
    ap.add_argument("--antigens", type=int, default=8)
    ap.add_argument("--groups", type=int, default=3)
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--beads", type=int, default=3)
    ap.add_argument("--reads", type=int, default=4000, help="reads per sample")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = Path(args.out_dir)
    reads_dir = out / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)

    library = build_library(rng, args.oligos, args.antigens)
    oligos = {r["peptide_id"]: r["oligo"] for r in library}
    with open(out / "peptide_table.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(library[0]))
        writer.writeheader()
        writer.writerows(library)

    with open(out / "adapters.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "read", "end", "sequence", "anchored"])
        writer.writerow(["leader_r1", "R1", "5", ADAPTER, "false"])

    # Boost one shared tile everywhere and an overlapping pair in the first group.
    shared = args.oligos // 2
    pair = (1, 2)

    samples, trid = [], 0
    for b in range(args.beads):
        trid += 1
        samples.append({"technical_replicate_id": trid,
                        "sample_name": f"beads_{b + 1}", "participant_ID": "",
                        "sample_role": "beads", "control_status": "beads_only",
                        "boost": []})
    for g in range(args.groups):
        for r in range(args.replicates):
            trid += 1
            boost = [shared] + (list(pair) if g == 0 else [])
            samples.append({"technical_replicate_id": trid,
                            "sample_name": f"D{g + 1}_{r + 1}",
                            "participant_ID": f"donor_{g + 1}",
                            "sample_role": "serum", "control_status": "empirical",
                            "boost": boost})

    rows = []
    for sample in samples:
        counts = {}
        remaining = args.reads
        for pid in sample["boost"]:
            counts[pid] = int(args.reads * 0.12)
            remaining -= counts[pid]
        for _ in range(remaining):
            pid = rng.randrange(args.oligos)
            counts[pid] = counts.get(pid, 0) + 1
        name = sample["sample_name"]
        r1 = reads_dir / f"{name}_R1.fastq.gz"
        write_reads(r1, oligos, counts, rng)
        rows.append({k: v for k, v in sample.items() if k != "boost"}
                    | {"fastq_filepath": str(r1.resolve())})

    with open(out / "sample_table.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"[fixture] {args.oligos} oligos, {len(rows)} samples in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
