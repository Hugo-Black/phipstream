#!/usr/bin/env nextflow
/*
Statistics and enrichment workflow.

Author: Jared G. Galloway
*/

// Enable Nextflow DSL 2.
nextflow.enable.dsl=2

// Import the edgeR and BEER scoring subworkflow.
include { edgeR_BEER_workflows } from './edgeR_BEER.nf'

/*
Automatically computed statistics with no required annotations.
*/

process counts_per_million {
    input: path phip_data
    output: path "cpm.phip"
    shell:
    """
    #!/usr/bin/env python3

    from phippery.normalize import counts_per_million
    from phippery.utils import *
    
    ds = load("$phip_data")
    counts_per_million(ds)
    dump(ds, "cpm.phip")
    """
}

process size_factors {
    input: path phip_data
    output: path "sf.phip"
    shell:
    """
    #!/usr/bin/env python3

    from phippery.normalize import size_factors
    from phippery.utils import *
    
    ds = load("$phip_data")
    size_factors(ds)
    dump(ds, "sf.phip") 
    """
}

/*
Optional statistics that require both annotations and enabled params flags.
*/

process cpm_fold_enrichment {
    input: path phip_data
    output: path "fold_enr.phip"
    when: params.run_cpm_enr_workflow
    shell:
    """
    #!/usr/bin/env python3

    from phippery.normalize import enrichment
    from phippery.utils import *
    
    ds = load("$phip_data")
    lib_ds = ds_query(ds, "control_status == 'library'")
    enrichment(ds, lib_ds, data_table="cpm")
    dump(ds, "fold_enr.phip") 
    """
}


process fit_predict_zscore {
    input: path phip_data
    output: path "fit_predict_zscore.phip"
    when: params.run_zscore_fit_predict || params.summarize_by_organism
    shell:
    """
    fit-predict-zscore.py \
        -ds ${phip_data} \
        -o fit_predict_zscore.phip
    """
}


process gamma_poisson_enrichment {
    input: path phip_data
    output: path "gamma_poisson.phip"
    when: params.run_gamma_poisson
    shell:
    """
    #!/usr/bin/env python3

    from phippery.modeling import gamma_poisson_model
    from phippery.utils import *

    ds = load("$phip_data")
    gamma_poisson_model(
        ds,
        starting_alpha=${params.gamma_poisson_starting_alpha},
        starting_beta=${params.gamma_poisson_starting_beta},
        trim_percentile=${params.gamma_poisson_trim_percentile},
        data_table="size_factors",
        inplace=True,
        new_table_name="gamma_poisson_mlxp",
    )
    dump(ds, "gamma_poisson.phip")
    """
}


/*
Merge phippery datasets through the xarray merge path.
*/

process merge_binary_datasets {    
    input:
    path all_phip_datasets
    output:
    path "merged.phip"
    shell:
    """
    phippery merge -o merged.phip '*.phip'
    """
}


workflow STATS {

    take: dataset
    main:

    // Compute baseline statistics that do not depend on annotations.
    dataset | \
        (counts_per_million & size_factors) | \
        mix | set { auto_stats_ch }

    if( params.run_edgeR | params.run_BEER )
        dataset | edgeR_BEER_workflows | set { edgeR_BEER_ch }
    else
        Channel.empty() | set { edgeR_BEER_ch }

    // Launch optional annotation-dependent statistics.
    cpm_fold_enrichment(counts_per_million.out) | set { cpm_fold_enr_ch }
    fit_predict_zscore(counts_per_million.out) | set { fit_pred_zscore_ch }
    gamma_poisson_enrichment(size_factors.out)  | set { gamma_poisson_ch }

    // Collect all statistic datasets and merge them.
    auto_stats_ch.concat(
        cpm_fold_enr_ch,
        fit_pred_zscore_ch,
        gamma_poisson_ch,
        edgeR_BEER_ch
    ) | collect | merge_binary_datasets

    emit:
    merge_binary_datasets.out
} 
