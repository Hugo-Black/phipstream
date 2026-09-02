#!/usr/bin/env bash
# Create a synthetic dataset and exercise the command line stages.
#
# The default path runs QC sheet creation, trimming, alignment, and
# prioritisation. Set BEER=1 to include the slower enrichment sampler.
#
# All files are written below a temporary directory that is removed on exit.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PHIPSTREAM_PYTHON:-python3}"
WORK="$(mktemp -d -t phipstream_test_XXXXXX)"
trap 'rm -rf "$WORK"' EXIT

fail() { echo "FAIL: $*" >&2 ; exit 1 ; }
rows() { $PYTHON -c "import sys; print(sum(1 for _ in open(sys.argv[1])) - 1)" "$1" ; }

echo "workdir $WORK"
"$PYTHON" "$ROOT/tests/make_fixture.py" --out-dir "$WORK/fixture"

"$ROOT/bin/phipstream" qc \
    --sample-table "$WORK/fixture/sample_table.csv" \
    --output       "$WORK/qc/samplesheet.csv"
[[ "$(rows "$WORK/qc/samplesheet.csv")" == "12" ]] \
    || fail "samplesheet should list 12 samples"

"$ROOT/bin/phipstream" trim \
    --sample-table         "$WORK/fixture/sample_table.csv" \
    --adapters             "$WORK/fixture/adapters.csv" \
    --out-dir              "$WORK/trimmed" \
    --trimmed-sample-table "$WORK/trimmed_table.csv" \
    --jobs 4
[[ -s "$WORK/trimmed_table.csv" ]] || fail "no trimmed sample table"

"$ROOT/bin/phipstream" align \
    --sample-table  "$WORK/trimmed_table.csv" \
    --peptide-table "$WORK/fixture/peptide_table.csv" \
    --out-dir       "$WORK/alignment" \
    --jobs 4
COUNTS="$WORK/alignment/counts/se_r1.local.csv"
[[ -s "$COUNTS" ]] || fail "no counts matrix"
"$PYTHON" - "$COUNTS" <<'PY'
import csv
import sys
rows = list(csv.reader(open(sys.argv[1])))
header, body = rows[0], rows[1:]
assert len(header) == 13, f"expected 12 sample columns, got {len(header) - 1}"
assert len(body) == 200, f"expected 200 oligos, got {len(body)}"
total = sum(int(v) for r in body for v in r[1:])
assert total > 40000, f"only {total} reads aligned, expected most of 48000"
print(f"[test] counts matrix {len(body)} by {len(header) - 1}, {total} reads placed")
PY

if [[ "${BEER:-0}" == "1" ]]; then
    "$ROOT/bin/phipstream" score \
        --counts        "$COUNTS" \
        --sample-table  "$WORK/trimmed_table.csv" \
        --peptide-table "$WORK/fixture/peptide_table.csv" \
        --out-dir       "$WORK/enrichment"
    POSTERIOR="$WORK/enrichment/beer_posterior.csv.gz"
    HITS="$WORK/enrichment/beer_hits.csv.gz"
else
    echo "[test] enrichment step skipped, set BEER=1 to include it"
    "$PYTHON" "$ROOT/tests/synthetic_calls.py" \
        --counts       "$COUNTS" \
        --sample-table "$WORK/trimmed_table.csv" \
        --out-dir      "$WORK/enrichment"
    POSTERIOR="$WORK/enrichment/beer_posterior.csv.gz"
    HITS="$WORK/enrichment/beer_hits.csv.gz"
fi

"$ROOT/bin/phipstream" prioritise \
    --posterior     "$POSTERIOR" \
    --hits          "$HITS" \
    --sample-table  "$WORK/trimmed_table.csv" \
    --peptide-table "$WORK/fixture/peptide_table.csv" \
    --out-dir       "$WORK/prioritised" \
    --role serum
"$PYTHON" - "$WORK/prioritised/prioritisation_summary.json" <<'PY'
import json
import sys
s = json.load(open(sys.argv[1]))
assert s["n_shortlist"] == 3, f"expected 3 shortlisted peptides, got {s['n_shortlist']}"
for name in ("groups", "replicates", "adjacent"):
    assert s["criteria"][name] > 0, f"criterion {name} never fired"
print(f"[test] shortlist {s['n_shortlist']}, criteria {s['criteria']}")
PY

# The scorers computed in Python, on the same counts the R path uses. true_gp is
# left out because it is the one that needs statsmodels.
for METHOD in arcsinh aitchison
do
    "$ROOT/bin/phipstream" score \
        --counts        "$COUNTS" \
        --sample-table  "$WORK/trimmed_table.csv" \
        --peptide-table "$WORK/fixture/peptide_table.csv" \
        --out-dir       "$WORK/score_$METHOD" \
        --method        "$METHOD" > /dev/null
    "$PYTHON" - "$WORK/score_$METHOD/enrichment_summary.json" "$METHOD" <<'CHECK'
import json
import sys
s = json.load(open(sys.argv[1]))
assert s["method"] == sys.argv[2], s["method"]
assert s["threshold"] == 3.5, f"expected the 3.5 z cutoff, got {s['threshold']}"
assert s["n_peptides"] == 200, s["n_peptides"]
print(f"[test] {sys.argv[2]}: {s['n_calls']} calls at z > {s['threshold']}")
CHECK
done

# The replicate rule can only lower a score, so it can never add calls.
"$ROOT/bin/phipstream" score \
    --counts        "$COUNTS" \
    --sample-table  "$WORK/trimmed_table.csv" \
    --peptide-table "$WORK/fixture/peptide_table.csv" \
    --out-dir       "$WORK/score_arcsinh_rep" \
    --method        arcsinh --replicate-rule > /dev/null
"$PYTHON" - "$WORK/score_arcsinh" "$WORK/score_arcsinh_rep" <<'CHECK'
import json
import sys
plain = json.load(open(sys.argv[1] + "/enrichment_summary.json"))
ruled = json.load(open(sys.argv[2] + "/enrichment_summary.json"))
assert ruled["replicate_rule"] is True, "replicate rule not recorded"
assert ruled["n_calls"] <= plain["n_calls"], (
    f"replicate rule raised calls, {plain['n_calls']} to {ruled['n_calls']}")
print(f"[test] replicate rule: {plain['n_calls']} calls down to {ruled['n_calls']}")
CHECK

# The binned generalized Poisson scorers, when statsmodels is available. The
# zero inflated variant can only lower a tail probability, so it can never call
# fewer peptides than the plain fit at the same threshold.
if "$PYTHON" -c "import statsmodels" 2>/dev/null
then
    for METHOD in larman_gp xu_zigp
    do
        "$ROOT/bin/phipstream" score \
            --counts        "$COUNTS" \
            --sample-table  "$WORK/trimmed_table.csv" \
            --peptide-table "$WORK/fixture/peptide_table.csv" \
            --out-dir       "$WORK/score_$METHOD" \
            --method        "$METHOD" --gp-bins 5 > /dev/null
    done
    "$PYTHON" - "$WORK/score_larman_gp" "$WORK/score_xu_zigp" <<'CHECK'
import json
import sys
plain = json.load(open(sys.argv[1] + "/enrichment_summary.json"))
zeroed = json.load(open(sys.argv[2] + "/enrichment_summary.json"))
for s in (plain, zeroed):
    assert s["threshold"] == 2.3, f"expected the 2.3 cutoff, got {s['threshold']}"
    assert s["gp_bins"] == 5, s["gp_bins"]
assert zeroed["n_calls"] >= plain["n_calls"], (
    f"zero inflation lowered calls, {plain['n_calls']} to {zeroed['n_calls']}")
print(f"[test] larman_gp {plain['n_calls']} calls, "
      f"xu_zigp {zeroed['n_calls']} at -log10 p > 2.3")
CHECK
else
    echo "[test] statsmodels absent, binned generalized Poisson scorers skipped"
fi

# A combined config drives both routes. Check the stage chain gains the FastQC
# stage and that rendering keeps parameter names case sensitive, since a
# lowercased run_BEER is accepted in silence and then ignored.
cat > "$WORK/combined.conf" <<CONF
[phipstream]
sample_table  = $WORK/fixture/sample_table.csv
peptide_table = $WORK/fixture/peptide_table.csv
adapters      = $WORK/fixture/adapters.csv
out_dir       = $WORK/combined
fastqc        = true

[workflow]
read_length = 151
run_edgeR   = true
run_BEER    = false
CONF

"$ROOT/bin/phipstream" stages "$WORK/combined.conf" --dry-run > "$WORK/plan.txt"
grep -q "run_fastqc.py" "$WORK/plan.txt" \
    || fail "fastqc = true should add the FastQC stage"
head -1 "$WORK/plan.txt" | grep -q "5 stages" \
    || fail "combined config should plan 5 stages"

"$ROOT/bin/phipstream" stages "$WORK/combined.conf" --workflow --dry-run > /dev/null
RENDERED="$WORK/combined/workflow.config"
[[ -f "$RENDERED" ]] || fail "--workflow should render $RENDERED"
grep -q "run_BEER = false" "$RENDERED" \
    || fail "rendered config lost the case of run_BEER"
grep -q "adapter_r1_5p_list = '" "$RENDERED" \
    || fail "adapters should carry over from the phipstream section"
grep -q "fastqc_dir = '$WORK/combined/fastqc'" "$RENDERED" \
    || fail "fastqc_dir should point at the FastQC stage output"
echo "[test] combined config plans 5 stages and renders $(grep -c '=' "$RENDERED") workflow params"

echo "PASS end_to_end.sh"
