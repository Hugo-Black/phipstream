#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s", type=str)
parser.add_argument("-o", type=str)
parser.add_argument("--run_zscore_fit_predict", type=str)
args = parser.parse_args()

sample_table = pd.read_csv(args.s, header=0)
if "fastq_filepath" not in sample_table.columns:
    raise KeyError("Must provide 'fastq_filepath' column in the table")

assert "sample_id" not in sample_table.columns, "Cannot include a column named 'sample_id' in the sample table"

sample_table["sample_id"] = list(range(len(sample_table)))
sample_table.set_index("sample_id", inplace=True)

# Add fastq_r2_filepath when input tables omit optional R2 paths.
if "fastq_r2_filepath" not in sample_table.columns:
    sample_table["fastq_r2_filepath"] = ""
sample_table["fastq_r2_filepath"] = sample_table["fastq_r2_filepath"].fillna("")

# Extra checks are required when Z-score fitting is enabled.
if args.run_zscore_fit_predict == "true":
    # The model needs control_status to identify control rows.
    msg = "Must provide a column 'control_status'"
    assert "control_status" in sample_table.columns.values, msg

    # Stop unless at least two beads-only controls are present.
    msg = "Must provide >1 samples labeled 'beads_only' under 'control_status'"
    assert (sample_table["control_status"] == "beads_only").sum() > 1, msg

sample_table.to_csv(args.o, index=True, na_rep="NA")

