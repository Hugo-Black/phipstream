#!/usr/bin/env nextflow

// Library sequence QC scans peptide oligos for GC extremes, homopolymers,
// dinucleotide entropy, tandem repeats, and amino acid repeat clusters. It writes
// feature files without consuming reads and runs only when params.run_library_qc
// is true.

nextflow.enable.dsl = 2

process scan_library_features {

    publishDir "${params.results}", mode: 'copy'

    input:
    path peptide_table

    output:
    path "library_qc/", emit: out_dir

    when:
    params.run_library_qc

    script:
    def proteinFlag = params.library_qc_protein_col ? "--peptide-protein-col ${params.library_qc_protein_col}" : ""
    def idFlag      = params.library_qc_id_col      ? "--peptide-id-col ${params.library_qc_id_col}"           : ""
    def orgFlag     = params.library_qc_organism_col? "--peptide-organism-col ${params.library_qc_organism_col}": ""
    def trim5Flag   = params.library_qc_trim5  > 0  ? "--trim5 ${params.library_qc_trim5}"                     : ""
    def trim3Flag   = params.library_qc_trim3  > 0  ? "--trim3 ${params.library_qc_trim3}"                     : ""
    def threshFlag  = params.library_qc_thresholds  ? "--thresholds '${params.library_qc_thresholds}'"          : ""
    """
    mkdir -p library_qc
    scan-library-features.py \\
        --peptide-table ${peptide_table} \\
        --peptide-seq-col ${params.library_qc_dna_col} \\
        ${proteinFlag} \\
        ${idFlag} \\
        ${orgFlag} \\
        ${trim5Flag} \\
        ${trim3Flag} \\
        ${threshFlag} \\
        --out-dir library_qc
    """
}

workflow LIBRARY_QC {
    take:
    peptide_table_ch

    main:
    scan_library_features(peptide_table_ch)
}
