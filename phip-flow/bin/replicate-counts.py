#!/usr/bin/env python
"""
Update counts for duplicated oligo sequences in a phippery dataset.

Peptide rows that share an oligo sequence receive the summed raw count across
that duplicate group. The dataset is modified in memory and then written by the
caller.
"""

import pandas as pd
import numpy as np
import phippery
from phippery.utils import load, dump, get_annotation_table
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()


def replicate_oligo_counts(ds, peptide_oligo_feature="Oligo"):
    """Sum raw counts across peptide rows that share an oligo sequence."""

    # Retained reference block for the older duplicate-oligo approach.
    # Locate oligo sequences that occur more than once.
    #pep_anno_table = get_annotation_table(ds, "peptide")
    #oligo_vc = pep_anno_table["Oligo"].value_counts()

    ## For each oligo sequence with multiple peptide rows.
    #for oligo, count in oligo_vc[oligo_vc > 1].items():
    #    replicate_idxs = pep_anno_table[
    #            pep_anno_table["Oligo"]==oligo
    #    ].index.values

    #    # Sum counts across duplicated oligo rows.
    #    rep_pep_sums = ds.counts.loc[replicate_idxs, :].sum(axis=0).values

    #    # Copy the summed counts back onto each duplicate row.
    #    ds.counts.loc[replicate_idxs, :] = np.tile(rep_pep_sums, (count, 1))

    # Load peptide annotations before grouping duplicated oligos.
    pep_anno_table = get_annotation_table(ds, "peptide")

    # Iterate through peptide groups that share one oligo sequence.
    for oligo_seq, pep_anno_table_oligo in pep_anno_table.groupby(peptide_oligo_feature):

        # Determine whether this oligo appears in more than one peptide row.
        if pep_anno_table_oligo.shape[0] == 1:

            # Leave unique oligos unchanged.
            continue

        # For duplicates, sum counts across the shared oligo rows.
        idxs = pep_anno_table_oligo.index.values
        # rep_pep_sums = ds.counts.loc[idxs, :].sum(axis=0).values

        # Assign that sum to every peptide row using the shared oligo.
        ds.counts.loc[idxs, :] = np.tile(
                ds.counts.loc[idxs, :].sum(axis=0).values,
                (pep_anno_table_oligo.shape[0], 1)
        )

ds = phippery.load(args.ds)
replicate_oligo_counts(ds, "oligo")
phippery.dump(ds, args.o)
