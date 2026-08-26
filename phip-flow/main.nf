/*
 * This source file is distributed under the GNU General Public License v3.0.
 */


/*
 * PhIP-Flow Nextflow workflow for common PhIP-Seq analysis tasks.
 *
 * Fred Hutchinson Cancer Research Center, Seattle WA.
 *
 * Jared Galloway
 * Kevin Sung
 * Sam Minot
 * Erick Matsen
 */

/*
 * Enable the Nextflow DSL 2 syntax.
 */
nextflow.enable.dsl = 2

// Dataset configs must provide the main input paths. Failing early here gives
// a clearer message than letting a missing default path break later.
params.sample_table  = null
params.peptide_table = null
params.reads_prefix  = "$launchDir"
params.results       = "$launchDir/results"

if (!params.sample_table || !params.peptide_table)
    error "set params.sample_table and params.peptide_table in a dataset config"

log.info """\
phipstream
================================
sample_table    : $params.sample_table
peptide_table   : $params.peptide_table
results         : $params.results

"""

/*
 * Import workflow modules.
 */
nextflow.enable.dsl=2

include { ALIGN } from './workflows/alignment.nf'
include { STATS } from './workflows/statistics.nf'
include { DSOUT } from './workflows/output.nf'
include { AGG }   from './workflows/aggregate.nf'

// Optional phipstream modules are enabled through their own params flags. They
// stay disabled by default to preserve upstream phip-flow behavior.
include { LIBRARY_QC } from './workflows/library_qc.nf'
include { CONTAMINATION_DIAGNOSTIC } from './workflows/contamination_diagnostic.nf'

workflow {

    // Run the core pipeline with explicit module calls. ALIGN emits multiple
    // channels, so naming the channel sent to STATS avoids ambiguous piping.
    ALIGN()
    STATS(ALIGN.out.final_output)
    DSOUT(STATS.out)
    AGG(DSOUT.out)

    // Library QC reads the peptide table directly and does not depend on ALIGN.
    if ( params.run_library_qc ) {
        LIBRARY_QC(Channel.fromPath(params.peptide_table))
    }

    // The contamination diagnostic rechecks external FastQC overrepresented
    // sequence output against the library and an optional source panel.
    if ( params.run_contamination_diagnostic ) {
        CONTAMINATION_DIAGNOSTIC(ALIGN.out.peptide_fasta)
    }
}
