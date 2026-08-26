#!/usr/bin/env python3

import os
from typing import List
import pandas as pd
import logging
from scipy.stats import gmean

# Method summary
#
# Start from epitope Z-scores for each replicate. Replicates belonging to the
# same sample are combined into an epitope-level table containing mean Z-score,
# hit state, discordance state, and whether the epitope is in the public list.
#
# Organism summaries keep the strongest non-overlapping epitopes for each sample.
# Overlap is detected with exact shared k-mers, so two peptides count as
# overlapping when they share any k-mer of the configured length.
#
# Final organism rows report hit and discordance counts plus maximum and mean EBS
# values, split across all epitopes and public epitopes.
#
# Inputs are staged by the surrounding Nextflow process according to the selected
# phip-flow module configuration.

class AggregatePhIP:

    def __init__(self):

        # Build the logger format and handler.
        self.logger = self.setup_logging()

        # Load the replicate-to-sample mapping.
        self.sample_mapping = self.read_sample_mapping()

        # Load the public epitope reference list.
        self.public_epitopes = self.read_public_epitopes()

        # Load peptide annotations used for grouping.
        self.peptide_mapping = self.read_peptide_mapping()

        # Read the maximum allowed peptide overlap.
        self.max_overlap = int("!{params.max_overlap}")
        self.logger.info(f"Maximum overlap: {self.max_overlap}")

        # Read the minimum Z-score threshold for hits.
        self.zscore_threshold = float("!{params.zscore_threshold}")
        self.logger.info(f"Z-score threshold: {self.zscore_threshold}")

        # Load replicate Z-scores.
        zscores_fp = "!{params.dataset_prefix}_zscore.csv.gz"
        self.logger.info(f"Reading in z-scores from: {zscores_fp}")
        assert os.path.exists(zscores_fp)
        self.zscores = pd.read_csv(zscores_fp, index_col=0)

        # Load edgeR hit calls when they were produced.
        edgeR_hits_fp = "!{params.dataset_prefix}_edgeR_hits.csv.gz"
        if os.path.exists(edgeR_hits_fp):
            self.logger.info(f"Reading in edgeR hits from: {edgeR_hits_fp}")
            self.edgeR_hits = pd.read_csv(
                edgeR_hits_fp,
                index_col=0,
                true_values=["TRUE", "True", "true"],
                false_values=["FALSE", "False", "false"],
                na_values=["NA", "N/A", "Na", "na", "n/a"]
            )
            self.has_edgeR_hits = True
        else:
            self.has_edgeR_hits = False

        # Collapse replicates into the current sample.
        self.logger.info("Grouping replicates by sample")
        self.sample_table = self.group_replicates()

        # Apply the overlap filter and record which peptides pass.
        self.sample_table = self.apply_max_overlap_filter()

        # Save the organism-level summary table.
        self.sample_table.to_csv("!{sample_id}.peptide.ebs.csv.gz", index=None)

        # Summarise filtered peptides by organism.
        self.logger.info("Grouping peptides by organism")
        self.organism_table = self.group_organisms()

        # Save to CSV
        self.logger.info("Writing organism-level outputs to CSV")
        self.organism_table.to_csv("!{sample_id}.organism.summary.csv.gz", index=None)

        self.logger.info("Done")

    def setup_logging(self) -> logging.Logger:
        """Create the task logger."""

        # Set up logging
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [aggregate_organisms] %(message)s'
        )
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        # Send log messages to standard output.
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)

        return logger

    def read_sample_mapping(self) -> pd.Series:
        """Load the column that maps replicate ids to sample ids."""

        # The sample annotation table supplies replicate grouping metadata.
        sample_mapping_fp = "!{params.dataset_prefix}_sample_annotation_table.csv.gz"
        self.logger.info(f"Reading in sample mapping from: {sample_mapping_fp}")
        assert os.path.exists(sample_mapping_fp)

        # Load the annotation table.
        df = pd.read_csv(sample_mapping_fp, index_col=0)
        self.logger.info(f"Sample mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # Use the configured grouping column when present.
        sample_grouping_col = "!{params.sample_grouping_col}"
        if len(sample_grouping_col) > 0:

            # Validate that the requested grouping column exists.
            msg = f"Column '{sample_grouping_col}' not found ({', '.join(df.columns.values)})"
            assert sample_grouping_col in df.columns.values, msg

            # Return the requested replicate-to-sample mapping.
            return df[sample_grouping_col]

        # Without grouping metadata, handle each replicate separately.
        else:

            # Map every replicate id to itself.
            return {
                int(replicate_id): str(replicate_id)
                for replicate_id in df.index.values
            }

    def read_peptide_mapping(self) -> pd.DataFrame:
        """Load peptide annotations and rename the configured grouping columns."""

        peptide_mapping_fp = "!{params.dataset_prefix}_peptide_annotation_table.csv.gz"
        self.logger.info(f"Reading in peptide mappings from: {peptide_mapping_fp}")
        assert os.path.exists(peptide_mapping_fp)

        # Read in the table
        df = pd.read_csv(peptide_mapping_fp, index_col=0)
        self.logger.info(f"Peptide mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # Map the user-provided names to controlled values
        mapping = {
            # Column used to group peptides by organism.
            "!{params.peptide_org_col}": "organism",
            # Column holding the protein sequence used for public epitope matching.
            "!{params.peptide_seq_col}": "seq"
        }

        # Check each configured source column.
        for cname in mapping.keys():

            # Confirm the peptide table contains this column.
            msg = f"Column '{cname}' not found ({', '.join(df.columns.values)})"
            assert cname in df.columns.values, msg

        # Rename configured columns to internal names.
        df = df.rename(columns=mapping)

        # Keep only fields required downstream.
        df = df.reindex(
            columns=list(mapping.values())
        )

        # Mark public epitopes after removing sequence after the first stop codon.
        df = df.assign(
            public=df["seq"].apply(
                lambda s: s.split("*")[0]
            ).isin(self.public_epitopes)
        )

        self.logger.info(f"Public Epitopes: {df['public'].sum():,} / {df.shape[0]:,}")

        # Add peptide length for summary logic.
        df = df.assign(
            peptide_length=lambda d: d["seq"].apply(len)
        )

        return df

    def read_public_epitopes(self) -> List[str]:
        """Read the public epitope sequence list."""

        # Load the public epitope table.
        df = pd.read_csv("!{public_epitopes_csv}")
        self.logger.info(f"Public epitope table has {df.shape[0]:,} rows")

        # Column containing public epitope sequences.
        public_epitopes_col = "peptide_translate"

        msg = f"Column not found: {public_epitopes_col} in ({', '.join(df.columns.values)})"
        assert public_epitopes_col in df.columns.values, msg

        # Remove sequence after a stop codon marker.
        return df[
            public_epitopes_col
        ].apply(
            lambda s: s.split("*")[0]
        ).tolist()

    def group_replicates(self) -> pd.DataFrame:
        """Collapse replicate Z-scores into the current sample table."""

        # Select replicates assigned to this sample.
        replicates = [
            rep_i
            for rep_i in self.zscores.columns.values
            if self.sample_mapping.get(int(rep_i)) == '!{sample_id}'
        ]

        self.logger.info(f"Filtering down to the {len(replicates):,} replicates for sample '!{sample_id}'")
        assert len(replicates) > 0

        # Keep only selected replicate columns.
        df = self.zscores.reindex(columns=replicates)
        # Match edgeR calls to the same replicate set.
        if self.has_edgeR_hits:
            self.edgeR_hits = self.edgeR_hits.reindex(columns=replicates)
            
        # Calculate sample-level epitope metrics.
        df = df.assign(
            n_replicates=len(replicates),
            EBS=df.mean(axis=1),
            hit=df.apply(self.classify_hit, axis=1),
            edgeR_hit=(
                self.edgeR_hits.apply(self.classify_edgeR_hit, axis=1)
                if self.has_edgeR_hits
                else None
            ),
            sample='!{sample_id}'
        ).reset_index(
        ).rename(
            columns=dict(index="peptide")
        ).drop(
            columns=replicates + (
                ["edgeR_hit"] if not self.has_edgeR_hits else []
            )
        )

        # Attach the public epitope flag.
        df = df.assign(
            public=df["peptide"].apply(int).apply(
                lambda i: self.peptide_mapping["public"][i]
            )
        )

        return df

    def classify_hit(self, r):
        """Classify a peptide as hit, non-hit, or discordant across replicates."""

        # Compare replicate Z-scores with the hit threshold.
        hit_vec = r > self.zscore_threshold

        # Convert replicate calls into a single state.
        if hit_vec.all():
            return "TRUE"
        elif not hit_vec.any():
            return "FALSE"
        else:
            return "DISCORDANT"

    def classify_edgeR_hit(self, r: pd.Series) -> str:
        """Classify edgeR replicate calls as hit, non-hit, or discordant."""

        # Remove missing edgeR values before classification.
        r = r.dropna()
        if len(r) == 0:  # all edgeR values were missing
            return "NA"

        # Determine the hit type
        if r.all():
            return "TRUE"
        elif not r.any():
            return "FALSE"
        else:
            return "DISCORDANT"

    def apply_max_overlap_filter(self) -> pd.DataFrame:
        """Run the overlap filter within each sample and organism group."""

        # Filter each sample and organism combination independently.
        df = pd.concat([
            self.apply_max_overlap_filter_sub(d)
            for _, d in self.sample_table.assign(
                organism=lambda d: d["peptide"].apply(
                    self.peptide_mapping["organism"].get
                )
            ).groupby(
                ["sample", "organism"]
            )
        ])

        return df

    def apply_max_overlap_filter_sub(
        self,
        df: pd.DataFrame
    ) -> pd.DataFrame:

        # Add peptide sequences for k-mer overlap checks.
        df = df.assign(
            seq=df["peptide"].apply(
                self.peptide_mapping["seq"].get
            ).apply(
                lambda s: s.rstrip("*")
            )
        )

        # Visit peptides from highest to lowest EBS.
        df = df.sort_values(by="EBS", ascending=False)

        # Track k-mers already claimed by accepted peptides.
        kmers_seen = set()

        # Store pass or fail decisions in row order.
        passes_filter = list()

        # Evaluate strongest binders first.
        for _, r in df.iterrows():

            # Build this peptide's k-mer set.
            row_kmers = set([
                r["seq"][n:(n + self.max_overlap)]
                for n in range(len(r["seq"]) - self.max_overlap)
            ])

            # Accept peptides whose k-mers have not appeared in accepted rows.
            passes_filter.append(len(row_kmers & kmers_seen) == 0)

            # If it passes
            if passes_filter[-1]:

                # Mark this peptide's k-mers as covered.
                kmers_seen |= row_kmers

        # Add the overlap filter decision to the table.
        df = df.assign(
            passes_filter=passes_filter
        )

        # Remove temporary sequence data.
        return (
            df
            .drop(columns=["seq"])
            .sort_index()
        )

    def group_organisms(self) -> pd.DataFrame:
        """Summarise sample-level peptide rows by organism."""

        # Build one summary per sample and organism combination.
        df = pd.concat([
            self.group_sample_organisms(d, sample, organism)
            for (sample, organism), d in self.sample_table.assign(
                organism=lambda d: d["peptide"].apply(
                    self.peptide_mapping["organism"].get
                )
            ).groupby(
                ["sample", "organism"]
            )
        ]).fillna(
            0
        )

        return df

    def group_sample_organisms(
        self,
        df: pd.DataFrame,
        sample: str,
        organism: str
    ) -> pd.DataFrame:

        """Summarise filtered peptides for one sample and one organism."""

        # Keep only peptides that passed the overlap filter.
        df = df.query("passes_filter")

        # Return hit counts and EBS summaries for all and public subsets.
        dat = pd.DataFrame([{
            "sample": sample,
            "organism": organism,
            **{
                k: v
                for label, d in [
                    ("all", df),
                    ("public", df.query("public")),
                    ("hits", df.query("hit == 'TRUE'")),
                ]
                if d.shape[0] > 0
                for k, v in [
                    (f"n_hits_{label}", (d["hit"] == "TRUE").sum()),
                    (f"n_discordant_{label}", (d["hit"] == "DISCORDANT").sum()),
                    (f"max_ebs_{label}", d["EBS"].max()),
                    (f"mean_ebs_{label}", d["EBS"].mean()),
                    (f"gmean_ebs_{label}", gmean(d["EBS"]))
                ] + (
                    [
                        (f"n_edgeR_hits_{label}", (d["edgeR_hit"] == "TRUE").sum()),
                        (f"n_edgeR_discordant_{label}", (d["edgeR_hit"] == "DISCORDANT").sum()),
                    ]
                    if self.has_edgeR_hits
                    else []
                )
                if k not in [
                    "n_hits_hits",
                    "n_discordant_hits",
                    "gmean_ebs_all",
                    "gmean_ebs_public"
                ]
            }
        }])

        return dat


AggregatePhIP()
