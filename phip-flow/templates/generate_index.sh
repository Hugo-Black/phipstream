#!/bin/bash

set -euo pipefail

FASTA=!{oligo_fasta}
CPUS=!{task.cpus}
ALIGNER=!{params.aligner}

mkdir peptide_index
case "$ALIGNER" in
    bowtie1)
        bowtie-build \
            --threads $CPUS \
            $FASTA \
            peptide_index/peptide
        ;;
    bowtie2)
        bowtie2-build \
            --threads $CPUS \
            $FASTA \
            peptide_index/peptide
        ;;
    *)
        echo "Unknown aligner: $ALIGNER (must be 'bowtie1' or 'bowtie2')" >&2
        exit 1
        ;;
esac

# Record the aligner so downstream steps can verify the index type.
echo "$ALIGNER" > peptide_index/aligner.txt
