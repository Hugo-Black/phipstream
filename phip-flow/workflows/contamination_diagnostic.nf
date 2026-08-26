#!/usr/bin/env nextflow

// Tier 0 contamination diagnostic for FastQC overrepresented sequences. The
// process compares unique sequences with the peptide library and optional
// contamination FASTA files, then writes unassigned sequences for optional BLAST.
// Enable it with params.run_contamination_diagnostic and point params.fastqc_dir
// at zipped or unpacked FastQC output. Later tier diagnostics are documented in
// wiki/plan/phip-flow-contamination-diagnostic-module.md but are not implemented
// in this workflow.

nextflow.enable.dsl = 2

process scan_fastqc_overrepresented_tier0 {

    publishDir "${params.results}/contamination_diagnostic", mode: 'copy'

    input:
    path peptide_fasta
    path fastqc_dir
    path contamination_sources_dir

    output:
    path "tier0/", emit: out_dir

    when:
    params.run_contamination_diagnostic

    script:
    def srcFlag = contamination_sources_dir.name == 'NO_SOURCES' \
        ? "" \
        : "--contamination-sources-dir ${contamination_sources_dir}"
    """
    mkdir -p tier0
    scan-fastqc-overrepresented.py \\
        --fastqc-dir ${fastqc_dir} \\
        --library-fasta ${peptide_fasta} \\
        ${srcFlag} \\
        --k ${params.contamination_kmer_size} \\
        --top-n ${params.fastqc_overrepresented_top_n} \\
        --out-dir tier0
    """
}

workflow CONTAMINATION_DIAGNOSTIC {
    take:
    peptide_fasta_ch

    main:
    if ( !params.fastqc_dir || !file(params.fastqc_dir).exists() ) {
        log.warn "CONTAMINATION_DIAGNOSTIC: params.fastqc_dir ('${params.fastqc_dir}') " +
                 "does not resolve; Tier 0 skipped. Set params.fastqc_dir to an existing " +
                 "FastQC output directory (containing *_fastqc.zip archives) to enable."
        return
    }

    def src_path = (params.contamination_sources_dir && file(params.contamination_sources_dir).exists())
        ? file(params.contamination_sources_dir)
        : file('NO_SOURCES')

    scan_fastqc_overrepresented_tier0(
        peptide_fasta_ch,
        Channel.value(file(params.fastqc_dir)),
        Channel.value(src_path)
    )
}
