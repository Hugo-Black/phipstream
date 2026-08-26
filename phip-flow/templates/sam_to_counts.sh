#!/bin/bash

set -euo pipefail
CPUS=!{task.cpus}

# Convert SAM to a sorted BAM file.
samtools view -u -@ ${CPUS} !{sam_file} | \
    samtools sort -@ ${CPUS} - > !{sample_id}.bam

# Re-sort into the final BAM filename.
samtools sort -@ ${CPUS} !{sample_id}.bam -o !{sample_id}.sorted 

# Replace the temporary BAM with the final sorted BAM.
mv !{sample_id}.sorted !{sample_id}.bam

# Build a BAM index.
samtools index -b !{sample_id}.bam

# Count per reference sequence. A paired library is tallied by fragment, taking
# the primary first mate of a proper pair, because both mates of a fragment land
# on the same oligo and counting segments would score it twice.
if [ "$(samtools view -c -f 1 !{sample_id}.bam)" -gt 0 ]; then
    samtools idxstats !{sample_id}.bam | grep -v '^\*' | cut -f1 > refs.txt
    samtools view -f 66 -F 2304 !{sample_id}.bam | cut -f3 > hits.txt
    awk 'NR==FNR{c[$1]++; next} {print $1"\t"(($1 in c)?c[$1]:0)}' \
        hits.txt refs.txt > !{sample_id}.counts
else
    samtools idxstats !{sample_id}.bam | \
        cut -f 1,3 | \
        sed "/^*/d" > !{sample_id}.counts
fi
