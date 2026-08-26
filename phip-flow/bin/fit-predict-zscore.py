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
args = parser.parse_args()

ds = load(args.ds)
beads_ds = ds_query(ds, "control_status == 'beads_only'")

zscore_ds = zscore(
    ds,
    beads_ds,
    data_table='cpm',
    min_Npeptides_per_bin=300,
    lower_quantile_limit=0.05,
    upper_quantile_limit=0.95,
    inplace=False,
    new_table_name='zscore'
)

dump(zscore_ds, args.o)
