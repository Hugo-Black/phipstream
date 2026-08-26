#!/usr/bin/env python3

from collections import defaultdict
import os
from typing import List
import pandas as pd
import logging


def setup_logging() -> logging.Logger:
    """Create the task logger."""

    # Build the logger format and handler.
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [split_samples] %(message)s'
    )
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Send log output to standard output.
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    return logger

logger = setup_logging()

# The sample annotation table supplies replicate grouping metadata.
sample_mapping_fp = "!{params.dataset_prefix}_sample_annotation_table.csv.gz"
logger.info(f"Reading in sample mapping from: {sample_mapping_fp}")
assert os.path.exists(sample_mapping_fp)

# Load the sample annotation table.
df = pd.read_csv(sample_mapping_fp, index_col=0)
logger.info(f"Sample mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

# Use the configured grouping column when present.
sample_grouping_col = "!{params.sample_grouping_col}"
if len(sample_grouping_col) > 0:

    # Validate that the grouping column exists.
    msg = f"Column '{sample_grouping_col}' not found ({', '.join(df.columns.values)})"
    assert sample_grouping_col in df.columns.values, msg

    # Write the unique sample names to disk.
    df.reindex(
        columns=[sample_grouping_col]
    ).drop_duplicates(
    ).to_csv(
        "sample_list",
        header=None,
        index=None
    )

# Fall back when no grouping column is configured.
else:

    # Treat each replicate id as its own sample.
    with open("sample_list", "w") as handle:
        handle.write(
            "\n".join(list(map(str, df.index.values)))
        )
