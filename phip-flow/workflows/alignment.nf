#!/usr/bin/env nextflow

// Enable Nextflow DSL 2.
nextflow.enable.dsl=2

/*
Check the sample table and write the normalized copy used downstream.
*/
process validate_sample_table {
    input: path samples
    output: path "validated_sample_table.csv"
    script:
    """
    validate-sample-table.py \
        -s $samples \
        -o validated_sample_table.csv \
        --run_zscore_fit_predict ${params.run_zscore_fit_predict}
    """  
}

/*
Check the peptide table and write the normalized copy used downstream.
*/
// The validator expects a peptide id column in the input table.
process validate_peptide_table{
    input: path peptides
    output: path "validated_peptide_table.csv"
    script:
    """
    validate-peptide-table.py \
        -p $peptides \
        -o validated_peptide_table.csv
    """
}

// Write peptide oligo metadata as a FASTA reference.
process generate_fasta_reference {
    input: path peptide_table
    output: path "peptides.fasta"
    script:
    """
    generate-fasta.py \
        -pt $peptide_table \
        -o peptides.fasta
    """
}


// Build the aligner index from the FASTA reference.
process generate_index {
    input:
    path "oligo_fasta"
    output:
    tuple val("peptide_ref"), path("peptide_index")
    shell:    
    template "generate_index.sh"
}


// Optional cutadapt process for raw R1 and optional R2 reads. Multiple adapter
// entries per end are passed as separate flags so mixed library geometries can
// be trimmed in a single process.
process trim_adapters {
    label 'mem_small'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz")

    when:
    params.run_adapter_trim

    script:
    def r1_5p_flags = params.adapter_r1_5p_list?.tokenize(',')?.collect { "-g '${it.trim()}'" }?.join(' ') ?: ''
    def r1_3p_flags = params.adapter_r1_3p_list?.tokenize(',')?.collect { "-a '${it.trim()}'" }?.join(' ') ?: ''
    def r2_5p_flags = params.adapter_r2_5p_list?.tokenize(',')?.collect { "-G '${it.trim()}'" }?.join(' ') ?: ''
    def r2_3p_flags = params.adapter_r2_3p_list?.tokenize(',')?.collect { "-A '${it.trim()}'" }?.join(' ') ?: ''
    def wildcards   = params.adapter_match_read_wildcards ? '--match-read-wildcards' : ''
    def min_len     = params.adapter_min_length
    def paired      = (r2.name != 'NO_R2')
    if ( paired ) {
        """
        cutadapt \\
            ${r1_5p_flags} ${r1_3p_flags} \\
            ${r2_5p_flags} ${r2_3p_flags} \\
            ${wildcards} \\
            --minimum-length ${min_len} \\
            -j ${task.cpus} \\
            -o trimmed_${sample_id}_R1.fastq.gz \\
            -p trimmed_${sample_id}_R2.fastq.gz \\
            ${r1} ${r2}
        """
    } else {
        // For single-end input, keep the tuple shape by emitting a NO_R2 file.
        """
        cutadapt \\
            ${r1_5p_flags} ${r1_3p_flags} \\
            ${wildcards} \\
            --minimum-length ${min_len} \\
            -j ${task.cpus} \\
            -o trimmed_${sample_id}_R1.fastq.gz \\
            ${r1}
        touch trimmed_${sample_id}_R2.fastq.gz   # placeholder so channel shape matches
        """
    }
}


// Align each sample against the reference index.
process short_read_alignment {
    label 'alignment_tool'
    input:
    tuple val(sample_id), path(index), path(respective_replicate_path), path(r2_path)
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    shell:
    template "short_read_alignment.sh"

}


// Generate alignment statistics for each SAM file.
process sam_to_stats {
    input:
    tuple val(sample_id), path(sam_file)
    output:
    path "${sample_id}.stats"
    shell:
    template "sam_to_stats.sh"
}


// Count alignments per peptide for each sample.
process sam_to_counts {
    input: tuple val(sample_id), path(sam_file)
    output: path "${sample_id}.counts"
    shell:
    template "sam_to_counts.sh"
}


// Merge count files, alignment stats, and metadata into one phippery dataset.
process collect_phip_data {
    input:
    path all_counts_files
    path all_alignment_stats
    path sample_table 
    path peptide_table 
    output:
    path "data.phip"

    shell:
    """
    merge-counts-stats.py \
        -st ${sample_table} \
        -pt ${peptide_table} \
        -cfp "*.counts" \
        -sfp "*.stats" \
        -o data.phip
    """
}

process replicate_counts {
    input: path ds
    output: path "replicated_counts.phip"
    script: 
    """
    replicate-counts.py \
        -ds ${ds} \
        -o replicated_counts.phip
    """
}

workflow ALIGN {

    main:
        sample_ch = Channel.fromPath(params.sample_table)
        peptide_ch = Channel.fromPath(params.peptide_table)

        validate_sample_table(sample_ch)
        validate_peptide_table(peptide_ch) \
            | generate_fasta_reference | generate_index

        // Use reads_prefix only for relative paths; absolute paths pass through.
        def read_file = { p ->
            p.startsWith('/') ? file(p, checkIfExists:true)
                              : file("${params.reads_prefix}/${p}", checkIfExists:true)
        }

        validate_sample_table.out
            .splitCsv(header:true )
            .map{ row ->
                def r2 = row.fastq_r2_filepath?.trim()
                def r2_file = (r2 && r2 != "" && r2 != "NA")
                    ? read_file(r2)
                    : file("NO_R2")
                tuple(
                    "peptide_ref",
                    row.sample_id,
                    read_file(row.fastq_filepath.trim()),
                    r2_file
                )
            }
            .set { samples_ch }

        // If trimming is enabled, send sample read tuples through cutadapt before
        // recreating the reference-tagged tuples consumed by alignment.
        if ( params.run_adapter_trim ) {
            samples_ch
                .map { ref, sample_id, r1, r2 -> tuple(sample_id, r1, r2) }
                .set { trim_in_ch }
            trim_adapters(trim_in_ch)
            samples_ch = trim_adapters.out.map { sample_id, r1, r2 ->
                tuple("peptide_ref", sample_id, r1, r2)
            }
        }

        short_read_alignment(
            generate_index.out
                .cross(samples_ch)
                .map{ ref, sample ->
                    tuple(
                        sample[1],          // sample id
                        file(ref[1]),       // index files
                        file(sample[2]),    // R1 path
                        sample[3],          // R2 path or NO_R2
                    )
                }
        ) | (sam_to_counts & sam_to_stats)

        ds = collect_phip_data(
            sam_to_counts.out.toSortedList(),
            sam_to_stats.out.toSortedList(),
            validate_sample_table.out,
            validate_peptide_table.out
        )

        final_output = ds
        if ( params.replicate_sequence_counts )
            final_output = replicate_counts(ds)

    emit:
        final_output
        peptide_fasta = generate_fasta_reference.out
}
