#!/usr/bin/env nextflow

// Enable Nextflow DSL 2.
nextflow.enable.dsl=2


// Build the sample list used to shard aggregate_organisms.
process split_samples {

    input:
        // Wide-format CSV outputs.
        path "*"

    output: path "sample_list"
    when: params.summarize_by_organism
    shell:
    template "split_samples.py"
}

process aggregate_organisms {
    tag "${sample_id}"
    cpus 1
    memory "4.GB"
    input:
        // Wide-format CSV outputs.
        tuple path("*"), val(sample_id)
        // Public epitope CSV input.
        path public_epitopes_csv
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template "aggregate_organisms.py"
}

process join_organisms {
    publishDir "$params.results/aggregated_data/", mode: 'copy', overwrite: true
    input: path "input/"
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template 'join_organisms.py'
}

workflow AGG {
    take:
        dump_binary
        dump_wide_csv
        dump_tall_csv
    main:

    // Create the list of samples to aggregate.
    split_samples(dump_wide_csv)

    aggregate_organisms(
        dump_wide_csv
            .toSortedList()
            .combine(
                split_samples
                    .out
                    .splitText(){it.replace("\n", "")}
            ),
        file("${params.public_epitopes_csv}")
    )

    join_organisms(
        aggregate_organisms
            .out
            .flatten()
            .toSortedList()
    )
}


