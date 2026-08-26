#!/usr/bin/env python3
"""Bowtie 2 alignment and count collection for trimmed reads.

Reads are aligned against nucleotide oligos from the peptide table. Mode selects
which reads are used and how the strand is constrained.

  se_r1  single-end on R1, restricted to the forward strand with --norc
  se_r2  single-end on R2, restricted to the reverse strand with --nofw,
         because R2 reads the reverse complement of the sense insert
  pe     paired-end, with the strand left to the fragment orientation

Preset `local` maps to bowtie2 --very-sensitive-local and allows clipped read
ends. Preset `end-to-end` maps to bowtie2 --very-sensitive and requires the full
read to align.

Paired-end counting tallies each fragment once, from the primary first mate of a
proper pair. Counting aligned segments instead would score both mates of every
fragment, and both mates land on the same oligo.

Outputs
  counts/<mode>.<preset>.csv   peptide_id rows by sample_name columns, zero filled
  align_summary.csv            per sample read totals and alignment rate
"""
import argparse
import csv
import gzip
import os
import shlex
import shutil
import statistics
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import containers  # noqa: E402

PRESETS = {"local": "--very-sensitive-local", "end-to-end": "--very-sensitive"}
# Strand constraint per mode. Paired-end takes none, because the pair itself
# fixes the orientation.
MODES = {"se_r1": "--norc", "se_r2": "--nofw", "pe": ""}
PAIRED = {"pe"}
USES_R2 = {"se_r2", "pe"}


def build_reference(peptide_table, out_dir):
    """Create the oligo FASTA and return its path with peptide order metadata."""
    ref_dir = Path(out_dir) / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    fasta = ref_dir / "peptide.fasta"
    ids, lengths = [], []
    with open(fasta, "w") as handle:
        for row in csv.DictReader(open(peptide_table)):
            pid = (row.get("peptide_id") or "").strip()
            oligo = (row.get("oligo") or "").strip().upper()
            if not pid or not oligo:
                continue
            handle.write(f">{pid}\n{oligo}\n")
            ids.append(pid)
            lengths.append(len(oligo))
    if not ids:
        sys.exit(f"no peptide_id/oligo rows found in {peptide_table}")
    median_len = int(statistics.median(lengths))
    return str(fasta), ids, median_len


def build_index(fasta, out_dir, threads, binds):
    index = str(Path(out_dir) / "reference" / "bt2" / "peptide")
    Path(index).parent.mkdir(parents=True, exist_ok=True)
    if Path(index + ".1.bt2").exists():
        return index
    print("[align] building bowtie2 index")
    cmd = (f"bowtie2-build --quiet --threads {threads} "
           f"{shlex.quote(fasta)} {shlex.quote(index)}")
    containers.run("bowtie2", cmd, binds, Path(index).parent)
    return index


def median_read_length(fastq, cap=50000):
    if not fastq or not os.path.exists(fastq):
        return ""
    lengths = []
    try:
        with gzip.open(fastq, "rt") as handle:
            for i, line in enumerate(handle):
                if i % 4 == 1:
                    lengths.append(len(line.rstrip("\n")))
                    if len(lengths) >= cap:
                        break
    except OSError:
        return ""
    return round(statistics.median(lengths), 1) if lengths else ""


def align_and_count(index, r1, r2, mode, preset, threads, scratch, binds):
    sam = os.path.join(scratch, "aligned.sam")
    bam = os.path.join(scratch, "aligned.bam")
    quoted = {k: shlex.quote(v) for k, v in
              {"index": index, "sam": sam, "bam": bam, "scratch": scratch}.items()}
    if mode == "pe":
        reads = f"-1 {shlex.quote(r1)} -2 {shlex.quote(r2)}"
    else:
        reads = f"-U {shlex.quote(r2 if mode == 'se_r2' else r1)}"
    strand = f"{MODES[mode]} " if MODES[mode] else ""
    align = (f"bowtie2 {PRESETS[preset]} {strand}-p {threads} "
             f"-x {quoted['index']} {reads} -S {quoted['sam']}")
    if containers.run("bowtie2", align, binds | {scratch}, scratch,
                      capture=True, check=False).returncode != 0:
        return None

    # -f 66 keeps the first mate of a proper pair, -F 2304 drops secondary and
    # supplementary records, so one fragment contributes one count.
    if mode in PAIRED:
        tally = (f"samtools idxstats {quoted['bam']} | grep -v '^\\*' | cut -f1 "
                 f"> {quoted['scratch']}/refs.txt && "
                 f"samtools view -f 66 -F 2304 {quoted['bam']} | cut -f3 "
                 f"> {quoted['scratch']}/hits.txt && "
                 f"awk 'NR==FNR{{c[$1]++; next}} {{print $1\"\\t\"(($1 in c)?c[$1]:0)}}' "
                 f"{quoted['scratch']}/hits.txt {quoted['scratch']}/refs.txt "
                 f"> {quoted['scratch']}/counts.txt")
    else:
        tally = (f"samtools idxstats {quoted['bam']} | cut -f1,3 | grep -v '^\\*' "
                 f"> {quoted['scratch']}/counts.txt")
    post = (f"samtools sort -@ {threads} -o {quoted['bam']} {quoted['sam']} && "
            f"samtools index {quoted['bam']} && "
            f"samtools flagstat {quoted['bam']} > {quoted['scratch']}/flagstat.txt && "
            f"{tally}")
    if containers.run("samtools", post, binds | {scratch}, scratch,
                      capture=True, check=False).returncode != 0:
        return None

    total = mapped = proper = 0
    for line in Path(scratch, "flagstat.txt").read_text().splitlines():
        tail = line.split(" + 0 ", 1)[-1]
        if tail.startswith("primary mapped"):
            mapped = int(line.split()[0])
        elif tail.startswith("properly paired"):
            proper = int(line.split()[0])
        elif tail.startswith("primary") and "duplicate" not in tail:
            total = int(line.split()[0])
    counts = {}
    for line in Path(scratch, "counts.txt").read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] and int(parts[1]):
            counts[parts[0]] = int(parts[1])
    return {"total": total, "mapped": mapped, "proper": proper, "counts": counts}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-table", required=True,
                    help="trimmed sample table with sample_name and fastq_filepath")
    ap.add_argument("--peptide-table", required=True,
                    help="peptide table with peptide_id and oligo columns")
    ap.add_argument("--out-dir", required=True,
                    help="destination for counts, the reference and the summary")
    ap.add_argument("--preset", choices=sorted(PRESETS), default="local",
                    help="local allows clipped read ends, end-to-end requires "
                         "the whole read to align")
    ap.add_argument("--mode", choices=sorted(MODES), default="se_r1",
                    help="reads to align. se_r1 uses R1 on the forward strand, "
                         "se_r2 uses R2 on the reverse strand, pe uses both "
                         "and counts each fragment once (default se_r1)")
    ap.add_argument("--jobs", type=int, default=4, help="samples aligned at once")
    ap.add_argument("--threads", type=int, default=2, help="threads per sample")
    ap.add_argument("--scratch", default="", help="directory for intermediate BAM files")
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta, oligo_ids, median_len = build_reference(args.peptide_table, out_dir)
    print(f"[align] reference: {len(oligo_ids)} oligos, median length {median_len}")

    rows = list(csv.DictReader(open(args.sample_table)))
    kept = [r for r in rows
            if (r.get("sample_role") or "").strip().lower() != "undetermined"]
    if len(kept) < len(rows):
        print(f"[align] skipped {len(rows) - len(kept)} undetermined row(s)")
    if not kept:
        sys.exit("empty sample table")
    if "sample_name" not in kept[0]:
        sys.exit("sample table needs a sample_name column")

    needs_r2 = args.mode in USES_R2
    if needs_r2 and "fastq_r2_filepath" not in kept[0]:
        sys.exit(f"mode {args.mode} needs a fastq_r2_filepath column in the "
                 "sample table")
    binds = {str(out_dir)}
    for row in kept:
        paths = [(row.get("fastq_filepath") or "").strip()]
        if needs_r2:
            r2 = (row.get("fastq_r2_filepath") or "").strip()
            if not r2:
                sys.exit(f"mode {args.mode} needs R2 for {row['sample_name']}")
            paths.append(r2)
        for value in paths:
            path = Path(value)
            if not path.exists():
                sys.exit(f"missing read file: {path}")
            binds.add(str(path.resolve().parent))
    scratch_root = args.scratch or None
    if scratch_root:
        Path(scratch_root).mkdir(parents=True, exist_ok=True)
        binds.add(str(Path(scratch_root).resolve()))

    index = build_index(fasta, out_dir, args.threads, binds)
    detail = " ".join(x for x in (PRESETS[args.preset], MODES[args.mode]) if x)
    print(f"[align] {len(kept)} samples, mode {args.mode}, preset "
          f"{args.preset} ({detail})")

    def job(row):
        r1 = str(Path((row.get("fastq_filepath") or "").strip()).resolve())
        r2v = (row.get("fastq_r2_filepath") or "").strip()
        r2 = str(Path(r2v).resolve()) if r2v else ""
        scratch = tempfile.mkdtemp(dir=scratch_root, prefix="phipstream_aln_")
        try:
            result = align_and_count(index, r1, r2, args.mode, args.preset,
                                     args.threads, scratch, binds)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        return row, (r2 if args.mode == "se_r2" else r1), result

    per_sample, summary, failures = {}, [], []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(job, r) for r in kept]
        for done, fut in enumerate(as_completed(futures), 1):
            row, r1, result = fut.result()
            name = row["sample_name"]
            if result is None:
                failures.append(name)
                continue
            per_sample[name] = result["counts"]
            total = result["total"]
            summary.append({
                "sample_name": name,
                "mode": args.mode,
                "preset": args.preset,
                "sample_role": (row.get("sample_role") or "").strip(),
                "total_reads": total,
                "mapped_reads": result["mapped"],
                "pct_aligned": 100.0 * result["mapped"] / total if total else "",
                "proper_pairs": result["proper"] if args.mode in PAIRED else "",
                "read_len_median": median_read_length(r1),
                "ref_len": median_len,
            })
            if done % 25 == 0:
                print(f"  [{done}/{len(kept)}] aligned")

    with open(out_dir / "align_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_name", "mode", "preset", "sample_role", "total_reads",
            "mapped_reads", "pct_aligned", "proper_pairs", "read_len_median",
            "ref_len"])
        writer.writeheader()
        writer.writerows(sorted(summary, key=lambda r: r["sample_name"]))

    counts_dir = out_dir / "counts"
    counts_dir.mkdir(exist_ok=True)
    columns = sorted(per_sample)
    counts_path = counts_dir / f"{args.mode}.{args.preset}.csv"
    seen, ordered_ids = set(), []
    for pid in oligo_ids:
        if pid not in seen:
            ordered_ids.append(pid)
            seen.add(pid)
    with open(counts_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peptide_id"] + columns)
        for pid in ordered_ids:
            writer.writerow([pid] + [per_sample[c].get(pid, 0) for c in columns])

    print(f"[align] wrote {counts_path} ({len(ordered_ids)} oligos "
          f"by {len(columns)} samples)")
    if failures:
        sys.exit(f"{len(failures)} sample(s) failed to align: {failures[:10]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
