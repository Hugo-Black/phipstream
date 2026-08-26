#!/usr/bin/env python3
"""Screen FastQC overrepresented sequences for likely contamination sources.

The script reads overrepresented-sequence tables produced by FastQC or
seqinspector and classifies the most common sequences with the same k-mer based
attribution used by the read-level contamination probe.

Accepted inputs are either unpacked `*_fastqc/fastqc_data.txt` directories or
FastQC `*_fastqc.zip` archives.

The classifier compares each sequence against the peptide library and optional
contamination-source FASTA files using 25-mers by default.

Outputs:
    overrepresented_pooled.csv       one row per unique sequence across samples
    overrepresented_per_sample.csv   per-sample top-N rows with classifications
    overrepresented_summary.json     bucket counts and representative sequences
    overrepresented_for_blast.fasta  unassigned sequences for optional BLAST

Usage:
    scan-fastqc-overrepresented.py \n        --fastqc-dir results/<dataset>_out/qc/fastqc \n        --library-fasta <path>/peptide_library.fasta \n        --contamination-sources-dir <path>/contamination_sources \n        --out-dir /tmp/tier0_demo
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from library_qc_lib import gc_content, parse_fasta  # noqa: E402


# ---------------------------------------------------------------------------
# FastQC table parsing
# ---------------------------------------------------------------------------

_OVERREP_HEADER_RE = re.compile(r"^>>Overrepresented sequences\b")
_END_MODULE_RE     = re.compile(r"^>>END_MODULE")


def _extract_overrep_from_text(text: str) -> list[dict]:
    """Extract overrepresented sequence rows from fastqc_data.txt content."""
    lines = text.splitlines()
    in_section = False
    rows: list[dict] = []
    for line in lines:
        if _OVERREP_HEADER_RE.search(line):
            in_section = True
            continue
        if in_section and _END_MODULE_RE.search(line):
            break
        if not in_section:
            continue
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        seq, count, pct, src = parts[0], parts[1], parts[2], "\t".join(parts[3:])
        try:
            rows.append({
                "sequence": seq.upper(),
                "count": int(count),
                "percentage": float(pct),
                "fastqc_source": src,
            })
        except ValueError:
            continue
    return rows


def _iter_fastqc_files(fastqc_dir: Path):
    """Yield sample identifiers and parsed rows for each FastQC result."""
    # First handle unpacked fastqc_data.txt outputs.
    for p in fastqc_dir.rglob("fastqc_data.txt"):
        sample = p.parent.name
        if sample.endswith("_fastqc"):
            sample = sample[: -len("_fastqc")]
        yield sample, _extract_overrep_from_text(p.read_text())

    # Then handle the standard zipped FastQC archives.
    for zp in fastqc_dir.rglob("*_fastqc.zip"):
        sample = zp.stem
        if sample.endswith("_fastqc"):
            sample = sample[: -len("_fastqc")]
        try:
            with zipfile.ZipFile(zp) as zf:
                members = [m for m in zf.namelist() if m.endswith("fastqc_data.txt")]
                if not members:
                    continue
                with zf.open(members[0]) as fh:
                    text = io.TextIOWrapper(fh, encoding="utf-8").read()
                yield sample, _extract_overrep_from_text(text)
        except zipfile.BadZipFile:
            continue


# ---------------------------------------------------------------------------
# K-mer based attribution
# ---------------------------------------------------------------------------

def build_kmer_set(seqs: list[str], k: int) -> set[str]:
    kmers: set[str] = set()
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - k + 1):
            kmers.add(s[i : i + k])
    return kmers


def kmer_classify(seq: str, lib: set[str], src_kmers: dict[str, set[str]], k: int) -> tuple[str, list[str]]:
    """Return a classification label and matched contamination source names.

    Labels are library, vector_<name>, both, or neither. If several sources
    match, the primary label uses the first match in lookup order and the full
    matched source list is returned separately.
    """
    s = seq.upper()
    if len(s) < k:
        return "neither", []
    has_lib = False
    matched_sources: list[str] = []
    for i in range(len(s) - k + 1):
        km = s[i : i + k]
        if km in lib:
            has_lib = True
        for name, src in src_kmers.items():
            if km in src and name not in matched_sources:
                matched_sources.append(name)
    if has_lib and matched_sources:
        return "both", matched_sources
    if has_lib:
        return "library", []
    if matched_sources:
        return f"vector_{matched_sources[0]}", matched_sources
    return "neither", []


# ---------------------------------------------------------------------------
# Contamination source loader
# ---------------------------------------------------------------------------

def _load_contamination_sources(sources_dir: Path) -> dict[str, list[str]]:
    """Load sequences from each FASTA, gzipped FASTA, or zip in a directory."""
    if not sources_dir or not sources_dir.exists():
        return {}
    out: dict[str, list[str]] = {}
    extensions = {".fasta", ".fa", ".fna"}
    for p in sorted(sources_dir.iterdir()):
        # Ignore hidden files and browser metadata sidecars.
        if p.name.startswith(".") or p.name.endswith("Zone.Identifier"):
            continue
        name = p.stem
        if p.suffix.lower() in extensions:
            seqs = list(parse_fasta(p).values())
            if seqs:
                out[name] = seqs
        elif p.suffix.lower() == ".gz" and p.stem.endswith((".fasta", ".fa", ".fna")):
            base = p.stem.split(".")[0]
            with gzip.open(p, "rt") as fh:
                seqs = list(_parse_fasta_text(fh.read()).values())
            if seqs:
                out[base] = seqs
        elif p.suffix.lower() == ".zip":
            # Read the first FASTA-like member found in the zip.
            try:
                with zipfile.ZipFile(p) as zf:
                    for member in zf.namelist():
                        if member.endswith(tuple(extensions)):
                            with zf.open(member) as fh:
                                text = io.TextIOWrapper(fh, encoding="utf-8").read()
                            seqs = list(_parse_fasta_text(text).values())
                            if seqs:
                                out[name] = seqs
                            break
            except zipfile.BadZipFile:
                continue
    return out


def _parse_fasta_text(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    header: str | None = None
    buf: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                out[header] = "".join(buf).upper()
            header = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if header is not None:
        out[header] = "".join(buf).upper()
    return out


# ---------------------------------------------------------------------------
# Command line entry point
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastqc-dir", type=Path, required=True,
                   help="Directory containing FastQC output (fastqc_data.txt or *_fastqc.zip)")
    p.add_argument("--library-fasta", type=Path, required=True,
                   help="Peptide library FASTA")
    p.add_argument("--contamination-sources-dir", type=Path, default=None,
                   help="Optional directory of contamination-source FASTAs / zips")
    p.add_argument("--k", type=int, default=25,
                   help="k-mer size (default 25)")
    p.add_argument("--top-n", type=int, default=50,
                   help="Top N overrepresented sequences per sample to retain")
    p.add_argument("--min-percentage", type=float, default=0.1,
                   help="Minimum FastQC percentage to retain (default 0.1)")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load k-mers from the peptide library.
    lib_seqs = list(parse_fasta(args.library_fasta).values())
    if not lib_seqs:
        raise SystemExit(f"ERROR: empty library FASTA at {args.library_fasta}")
    lib_kmers = build_kmer_set(lib_seqs, args.k)

    # Load k-mers from optional contamination source files.
    contam = _load_contamination_sources(args.contamination_sources_dir) if args.contamination_sources_dir else {}
    src_kmers = {name: build_kmer_set(seqs, args.k) for name, seqs in contam.items()}

    print(f"Library k-mers: {len(lib_kmers):,} ({args.k}-mers)", file=sys.stderr)
    for name, kms in src_kmers.items():
        print(f"  source {name}: {len(kms):,} k-mers", file=sys.stderr)

    # Parse FastQC outputs and classify retained sequences.
    per_sample_rows: list[dict] = []
    seq_pool: dict[str, dict] = {}  # aggregate records keyed by sequence

    n_samples = 0
    n_total_overrep = 0
    for sample, rows in _iter_fastqc_files(args.fastqc_dir):
        if not rows:
            continue
        n_samples += 1
        rows = [r for r in rows if r["percentage"] >= args.min_percentage]
        rows = sorted(rows, key=lambda r: -r["percentage"])[: args.top_n]
        for r in rows:
            seq = r["sequence"]
            klass, src_match = kmer_classify(seq, lib_kmers, src_kmers, args.k)
            per_sample_rows.append({
                "sample_id": sample,
                "sequence": seq,
                "length": len(seq),
                "count_in_sample": r["count"],
                "percentage_in_sample": r["percentage"],
                "fastqc_source": r["fastqc_source"],
                "kmer_classification": klass,
                "kmer_sources_matched": ",".join(src_match),
                "gc_pct": round(100.0 * gc_content(seq), 2),
            })
            n_total_overrep += 1

            agg = seq_pool.setdefault(seq, {
                "sequence": seq,
                "length": len(seq),
                "gc_pct": round(100.0 * gc_content(seq), 2),
                "n_samples_with_hit": 0,
                "total_count": 0,
                "max_percentage": 0.0,
                "kmer_classification": klass,
                "kmer_sources_matched": ",".join(src_match),
                "fastqc_source_modal": Counter(),
            })
            agg["n_samples_with_hit"] += 1
            agg["total_count"] += r["count"]
            agg["max_percentage"] = max(agg["max_percentage"], r["percentage"])
            agg["fastqc_source_modal"][r["fastqc_source"]] += 1

    if n_samples == 0:
        raise SystemExit(
            f"ERROR: no FastQC output found under {args.fastqc_dir}. "
            "Expected *_fastqc.zip archives or unpacked fastqc_data.txt files."
        )

    # Convert aggregate records into output rows.
    pool_rows: list[dict] = []
    for seq, agg in seq_pool.items():
        modal = agg["fastqc_source_modal"].most_common(1)[0][0]
        pool_rows.append({
            "sequence": seq,
            "length": agg["length"],
            "gc_pct": agg["gc_pct"],
            "n_samples_with_hit": agg["n_samples_with_hit"],
            "total_count": agg["total_count"],
            "max_percentage": round(agg["max_percentage"], 3),
            "kmer_classification": agg["kmer_classification"],
            "kmer_sources_matched": agg["kmer_sources_matched"],
            "fastqc_source_modal": modal,
        })
    pool_rows.sort(key=lambda r: (-r["n_samples_with_hit"], -r["total_count"]))

    # Write the pooled and per-sample CSV outputs.
    pool_path = args.out_dir / "overrepresented_pooled.csv"
    with open(pool_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pool_rows[0].keys()))
        w.writeheader()
        w.writerows(pool_rows)

    sample_path = args.out_dir / "overrepresented_per_sample.csv"
    with open(sample_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_sample_rows[0].keys()))
        w.writeheader()
        w.writerows(per_sample_rows)

    # Write unassigned sequences for optional BLAST follow-up.
    blast_fa = args.out_dir / "overrepresented_for_blast.fasta"
    n_blast_candidates = 0
    with open(blast_fa, "w") as fh:
        for i, row in enumerate(pool_rows):
            if row["kmer_classification"] == "neither":
                fh.write(f">over_{i:05d} n_samples={row['n_samples_with_hit']} total={row['total_count']}\n")
                fh.write(f"{row['sequence']}\n")
                n_blast_candidates += 1

    # Write a compact JSON summary.
    bucket_counts: Counter[str] = Counter(r["kmer_classification"] for r in pool_rows)
    top_per_bucket: dict[str, dict] = {}
    for row in pool_rows:
        b = row["kmer_classification"]
        if b not in top_per_bucket:
            top_per_bucket[b] = {
                "sequence": row["sequence"],
                "n_samples_with_hit": row["n_samples_with_hit"],
                "total_count": row["total_count"],
                "fastqc_source_modal": row["fastqc_source_modal"],
            }
    summary = {
        "n_samples_with_overrep": n_samples,
        "n_total_overrep_rows": n_total_overrep,
        "n_unique_sequences": len(pool_rows),
        "kmer_classification_bucket_counts": dict(bucket_counts),
        "top_sequence_per_bucket": top_per_bucket,
        "n_blast_candidates": n_blast_candidates,
        "library_fasta": str(args.library_fasta.resolve()),
        "contamination_sources_dir": str(args.contamination_sources_dir.resolve()) if args.contamination_sources_dir else None,
        "contamination_sources_loaded": list(src_kmers.keys()),
        "k": args.k,
    }
    sum_path = args.out_dir / "overrepresented_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2))

    # Print a short run summary.
    print()
    print(f"FastQC samples processed: {n_samples}")
    print(f"Total overrepresented rows: {n_total_overrep}")
    print(f"Unique pooled sequences: {len(pool_rows)}")
    print(f"Bucket counts: {dict(bucket_counts)}")
    print(f"BLAST candidates (neither bucket): {n_blast_candidates}")
    print(f"Pooled CSV: {pool_path}")
    print(f"Per-sample CSV: {sample_path}")
    print(f"Summary JSON: {sum_path}")
    print(f"BLAST FASTA: {blast_fa}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
