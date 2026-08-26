"""Sequence feature utilities for the standalone library QC module.

The implementation uses only the Python standard library, allowing it to run in
the slim phip-flow container without adding Biopython.

Used by phip-flow/bin/scan-library-features.py.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# FASTA parsing and sequence loading
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> dict[str, str]:
    """Load FASTA records into a dict keyed by header."""
    out: dict[str, str] = {}
    header: str | None = None
    buf: list[str] = []
    with open(path) as fh:
        for raw in fh:
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
# DNA sequence metrics
# ---------------------------------------------------------------------------

def gc_content(seq: str) -> float:
    """Return the GC fraction for the whole sequence."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def gc_content_window(seq: str, w: int = 30) -> list[float]:
    """Return GC fractions for each window of width ``w``."""
    if len(seq) < w:
        return []
    out: list[float] = []
    # Seed the count from the first window.
    gc = sum(1 for b in seq[:w] if b in "GCgc")
    out.append(gc / w)
    # Slide the window while updating the GC count in place.
    for i in range(w, len(seq)):
        if seq[i - w] in "GCgc":
            gc -= 1
        if seq[i] in "GCgc":
            gc += 1
        out.append(gc / w)
    return out


def homopolymer_runs(seq: str) -> dict[str, int]:
    """Return the longest homopolymer run for each DNA base."""
    longest = {"A": 0, "C": 0, "G": 0, "T": 0}
    if not seq:
        return longest
    cur_base = seq[0].upper()
    cur_len = 1
    for b in seq[1:]:
        bu = b.upper()
        if bu == cur_base:
            cur_len += 1
        else:
            if cur_base in longest and cur_len > longest[cur_base]:
                longest[cur_base] = cur_len
            cur_base = bu
            cur_len = 1
    if cur_base in longest and cur_len > longest[cur_base]:
        longest[cur_base] = cur_len
    return longest


def dinucleotide_entropy(seq: str, w: int = 32) -> float:
    """Return the lowest dinucleotide Shannon entropy seen in any window.

    Low entropy marks repetitive sequence. The 32 bp window and 2.5 bit floor
    match the validation convention used for Twist handoff QC.
    """
    if len(seq) < w + 1:
        return float("inf")
    min_h = float("inf")
    for i in range(0, len(seq) - w):
        win = seq[i : i + w]
        counts: Counter[str] = Counter()
        for j in range(len(win) - 1):
            counts[win[j : j + 2]] += 1
        total = sum(counts.values())
        if total == 0:
            continue
        h = 0.0
        for c in counts.values():
            p = c / total
            h -= p * math.log2(p)
        if h < min_h:
            min_h = h
    return min_h


def tandem_repeat_span(seq: str, max_period: int = 12, min_copies: int = 3) -> int:
    """Return the base-pair span of the longest tandem repeat.

    Each period through max_period is tested at each start position, then
    extended until the repeated unit breaks. Only repeats with at least
    ``min_copies`` copies are counted.
    """
    n = len(seq)
    best = 0
    for period in range(1, max_period + 1):
        i = 0
        while i + period * min_copies <= n:
            unit = seq[i : i + period]
            j = i + period
            while j + period <= n and seq[j : j + period] == unit:
                j += period
            span = j - i
            copies = span // period
            if copies >= min_copies and span > best:
                best = span
            if j == i + period:
                i += 1
            else:
                i = j  # jump beyond the repeat that was just counted
    return best


# ---------------------------------------------------------------------------
# Protein sequence metrics
# ---------------------------------------------------------------------------

# Canonical codon table using NCBI translation table 1.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(dna: str, frame: int = 0) -> str:
    """Translate DNA in frame 0, 1, or 2; unknown codons become 'X'."""
    out: list[str] = []
    for i in range(frame, len(dna) - 2, 3):
        codon = dna[i : i + 3].upper()
        if len(codon) < 3:
            break
        out.append(_CODON_TABLE.get(codon, "X"))
    return "".join(out)


def aa_repeat_clusters(protein: str, min_len: int = 6) -> list[dict]:
    """Find mono-residue and di-residue low-complexity clusters.

    A cluster is a maximal contiguous span with no more than two amino acid
    symbols and at least 90 percent of residues coming from that small alphabet.
    Returned records include start, end, length, dominant_aa, motif_period, and
    the matched substring.

    These regions catch glycine-rich, alanine-rich, GA-repeat, and related
    sequences that are likely to be difficult during oligo synthesis.
    """
    if not protein:
        return []
    clusters: list[dict] = []
    n = len(protein)
    # Record maximal runs made from one amino acid.
    i = 0
    while i < n:
        j = i
        while j < n and protein[j] == protein[i] and protein[j] != "*":
            j += 1
        if j - i >= min_len:
            clusters.append({
                "start": i,
                "end": j,
                "length": j - i,
                "dominant_aa": protein[i],
                "motif_period": 1,
                "substring": protein[i:j],
            })
        i = j if j > i else i + 1

    # Search fixed windows for two-residue alphabets, then extend each match.
    i = 0
    while i < n - min_len:
        win = protein[i : i + min_len]
        distinct = set(win) - {"*"}
        if len(distinct) == 2:
            # Extend while every residue remains inside the two-letter alphabet.
            ab = distinct
            j = i + min_len
            while j < n and protein[j] in ab and protein[j] != "*":
                j += 1
            length = j - i
            # Keep clusters that are not already covered by a single-residue run.
            substr = protein[i:j]
            dom_a, dom_b = sorted(ab, key=lambda x: substr.count(x), reverse=True)
            if length >= min_len and any(
                (c["start"] <= i and c["end"] >= j) for c in clusters
            ):
                # A mono-residue cluster already accounts for this span.
                i += 1
                continue
            if length >= min_len:
                clusters.append({
                    "start": i,
                    "end": j,
                    "length": length,
                    "dominant_aa": f"{dom_a}{dom_b}",
                    "motif_period": 2,
                    "substring": substr,
                })
            i = j
        else:
            i += 1
    return clusters


# ---------------------------------------------------------------------------
# Per-oligo features and predicted synthesis failure
# ---------------------------------------------------------------------------

DEFAULT_THRESHOLDS = {
    "gc_window_max_extreme":           0.85,   # flag any window above this value
    "gc_window_max_high":              0.75,   # paired with aa_repeat_span
    "aa_repeat_span_min":              9,      # minimum low-complexity amino acid span
    "homopolymer_max":                 12,     # longest allowed single-base run
    "dinucleotide_entropy_min":        2.5,    # bit floor paired with tandem repeat span
    "tandem_repeat_span_min":          30,     # base-pair span floor
}


def featurize_oligo(
    dna: str,
    protein: str | None = None,
    *,
    trim5: int = 0,
    trim3: int = 0,
    gc_window: int = 30,
) -> dict:
    """Compute all sequence features for one oligo.

    ``trim5`` and ``trim3`` remove flanking adapter bases before DNA metrics are
    calculated. Protein metrics use the supplied protein sequence when present,
    otherwise the trimmed DNA is translated locally.
    """
    if trim5 < 0 or trim3 < 0:
        raise ValueError("trim5/trim3 must be non-negative")
    insert = dna[trim5 : len(dna) - trim3 if trim3 > 0 else len(dna)]

    gc_overall = gc_content(insert)
    gc_w = gc_content_window(insert, gc_window)
    if gc_w:
        gc_window_max = max(gc_w)
        gc_window_min = min(gc_w)
        gc_window_mean = sum(gc_w) / len(gc_w)
    else:
        gc_window_max = gc_window_min = gc_window_mean = gc_overall

    homo = homopolymer_runs(insert)
    dn_h = dinucleotide_entropy(insert, w=32) if len(insert) >= 33 else float("inf")
    tr_span = tandem_repeat_span(insert)

    if protein is None or not str(protein).strip():
        protein = translate(insert)
    # Keep only the protein sequence before the first stop codon.
    if "*" in protein:
        protein = protein.split("*")[0]

    clusters = aa_repeat_clusters(protein)
    if clusters:
        biggest = max(clusters, key=lambda c: c["length"])
        aa_span = biggest["length"]
        aa_class = biggest["dominant_aa"]
    else:
        aa_span = 0
        aa_class = ""

    return {
        "dna_len": len(insert),
        "protein_len": len(protein),
        "gc_overall": gc_overall,
        "gc_window_max": gc_window_max,
        "gc_window_min": gc_window_min,
        "gc_window_mean": gc_window_mean,
        "homopolymer_max_A": homo["A"],
        "homopolymer_max_C": homo["C"],
        "homopolymer_max_G": homo["G"],
        "homopolymer_max_T": homo["T"],
        "homopolymer_max": max(homo.values()),
        "dinucleotide_entropy_min": dn_h,
        "tandem_repeat_span_dna": tr_span,
        "aa_repeat_span": aa_span,
        "aa_repeat_class": aa_class,
        "protein": protein,
    }


def predicted_fail(features: dict, thresholds: dict | None = None) -> tuple[bool, list[str]]:
    """Return the predicted synthesis failure flag and its reasons.

    Reasons accumulate across thresholds. An empty reason list means the oligo
    passes this screen.
    """
    t = dict(DEFAULT_THRESHOLDS)
    if thresholds:
        t.update(thresholds)

    reasons: list[str] = []

    if features["gc_window_max"] >= t["gc_window_max_extreme"]:
        reasons.append("extreme_gc")
    if (
        features["gc_window_max"] >= t["gc_window_max_high"]
        and features["aa_repeat_span"] >= t["aa_repeat_span_min"]
    ):
        reasons.append("high_gc_low_complexity_aa_repeat")
    if features["homopolymer_max"] >= t["homopolymer_max"]:
        reasons.append("long_homopolymer")
    if (
        features["dinucleotide_entropy_min"] < t["dinucleotide_entropy_min"]
        and features["tandem_repeat_span_dna"] >= t["tandem_repeat_span_min"]
    ):
        reasons.append("tandem_repeat")

    return (bool(reasons), reasons)
