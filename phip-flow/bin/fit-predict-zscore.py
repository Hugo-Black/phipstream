#!/usr/bin/env python
"""Fit the phippery Z-score model and write the predicted dataset.

The input dataset must include sample annotations with control_status values.
Rows marked beads_only are used for fitting, and empirical rows are predicted
after peptide-level model fitting. The dataset must also contain the cpm table
created during the statistics step.

Background on the method is available in the phippery documentation:
https://matsengrp.github.io/phippery/
"""

from phippery.utils import * 
from phippery.modeling import zscore

import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
parser.add_argument(
    "--min-peptides-per-bin", type=int, default=300,
    help="smallest number of peptides a bin may hold. Peptides are ordered by "
         "their summed abundance across the beads-only samples and merged into "
         "bins until each reaches this floor, so the bin count is about the "
         "library size divided by this value. The default of 300 suits a "
         "library of order 100,000 peptides. A small library reaches the floor "
         "in a handful of bins, at which point a peptide is measured against "
         "others of quite different abundance (default 300)")
args = parser.parse_args()

ds = load(args.ds)
beads_ds = ds_query(ds, "control_status == 'beads_only'")

zscore_ds = zscore(
    ds,
    beads_ds,
    data_table='cpm',
    min_Npeptides_per_bin=args.min_peptides_per_bin,
    lower_quantile_limit=0.05,
    upper_quantile_limit=0.95,
    inplace=False,
    new_table_name='zscore'
)

dump(zscore_ds, args.o)
