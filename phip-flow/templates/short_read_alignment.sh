#!/bin/bash

: '
Align reads to the peptide library after positional trimming to the configured
library tile length.

Supported aligners come from params.aligner:
    bowtie1  upstream default with seed-length matching
    bowtie2  configured through params.bowtie2_optional_args

Alignment mode comes from params.alignment_mode:
    auto   use PE when R2 is present, otherwise SE_R1
    PE     require paired-end input
    SE_R1  single-end on R1, forward strand only
    SE_R2  single-end on R2, reverse strand only, since R2 reads the reverse
           complement of the sense insert
    SE     accepted as a name for SE_R1

The strand flag is added here rather than in the aligner argument strings,
because it depends on which read is being aligned.
'

set -euo pipefail

STREAM_FILE_CMD=!{params.fastq_stream_func}
FASTQ_R1=!{respective_replicate_path}
FASTQ_R2=!{r2_path}
INDEX=!{index}/peptide
ALIGN_OUT_FN=!{sample_id}.sam
READ_LENGTH=!{params.read_length}
PEPTIDE_LENGTH=!{params.oligo_tile_length}
CPUS=!{task.cpus}
MM=!{params.n_mismatches}
ALIGNER=!{params.aligner}
TRIM5=!{params.trim5}
ALIGNMENT_MODE=!{params.alignment_mode}

# 5 prime trimming shortens the read span available for alignment.
EFFECTIVE_READ_LENGTH=$(( READ_LENGTH - TRIM5 ))
if [ ${PEPTIDE_LENGTH} -lt ${EFFECTIVE_READ_LENGTH} ]; then
    let TRIM3=${EFFECTIVE_READ_LENGTH}-${PEPTIDE_LENGTH}
else
    TRIM3=0
fi

# Choose paired-end or single-end alignment mode.
R2_PRESENT=0
if [ -f "$FASTQ_R2" ] && [ "$FASTQ_R2" != "NO_R2" ]; then
    R2_PRESENT=1
fi

case "$ALIGNMENT_MODE" in
    auto)
        if [ $R2_PRESENT -eq 1 ]; then MODE=PE; else MODE=SE_R1; fi
        ;;
    PE|SE_R2)
        if [ $R2_PRESENT -ne 1 ]; then
            echo "alignment_mode=$ALIGNMENT_MODE but R2 is missing for !{sample_id}" >&2
            exit 1
        fi
        MODE=$ALIGNMENT_MODE
        ;;
    SE|SE_R1)
        MODE=SE_R1
        ;;
    *)
        echo "Unknown alignment_mode: $ALIGNMENT_MODE (auto, PE, SE_R1 or SE_R2)" >&2
        exit 1
        ;;
esac

# Read to align and strand constraint for single-end modes.
if [ "$MODE" = "SE_R2" ]; then
    SE_READ=$FASTQ_R2
    BT1_STRAND="--nofw"
    BT2_STRAND="--nofw"
else
    SE_READ=$FASTQ_R1
    BT1_STRAND="--norc"
    BT2_STRAND="--norc"
fi

case "$ALIGNER" in
    bowtie1)
        OP_ARGS="!{params.bowtie_optional_args}"
        echo "$OP_ARGS  (mode=$MODE)"
        if [ "$MODE" = "PE" ]; then
            echo "Running paired-end bowtie 1 for !{sample_id}"
            bowtie \
              --trim5 $TRIM5 \
              --trim3 $TRIM3 \
              --threads $CPUS \
              -n $MM \
              -l $PEPTIDE_LENGTH \
              $OP_ARGS \
              -x $INDEX \
              -1 <($STREAM_FILE_CMD $FASTQ_R1) \
              -2 <($STREAM_FILE_CMD $FASTQ_R2) \
              > $ALIGN_OUT_FN
        else
            echo "Running single-end bowtie 1 for !{sample_id}"
            $STREAM_FILE_CMD $SE_READ | bowtie \
              --trim5 $TRIM5 \
              --trim3 $TRIM3 \
              --threads $CPUS \
              -n $MM \
              -l $PEPTIDE_LENGTH \
              $OP_ARGS $BT1_STRAND \
              -x $INDEX - > $ALIGN_OUT_FN
        fi
        ;;
    bowtie2)
        OP_ARGS="!{params.bowtie2_optional_args}"
        echo "$OP_ARGS  (mode=$MODE)"
        if [ "$MODE" = "PE" ]; then
            echo "Running paired-end bowtie 2 for !{sample_id}"
            bowtie2 \
              --trim5 $TRIM5 \
              --trim3 $TRIM3 \
              --threads $CPUS \
              $OP_ARGS \
              -x $INDEX \
              -1 <($STREAM_FILE_CMD $FASTQ_R1) \
              -2 <($STREAM_FILE_CMD $FASTQ_R2) \
              -S $ALIGN_OUT_FN
        else
            echo "Running single-end bowtie 2 for !{sample_id}"
            bowtie2 \
              --trim5 $TRIM5 \
              --trim3 $TRIM3 \
              --threads $CPUS \
              $OP_ARGS $BT2_STRAND \
              -x $INDEX \
              -U <($STREAM_FILE_CMD $SE_READ) \
              -S $ALIGN_OUT_FN
        fi
        ;;
    *)
        echo "Unknown aligner: $ALIGNER (must be 'bowtie1' or 'bowtie2')" >&2
        exit 1
        ;;
esac
