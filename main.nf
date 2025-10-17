#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    jacgonisa/charla_nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/jacgonisa/charla_nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHARLA_NF } from './workflows/charla_nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_charla_nf_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_charla_nf_pipeline'

/*


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //

    // Backward compatibility: convert old col/ler params to new parent params
    if (params.col_bed_file && !params.parent1_bed) {
        params.parent1_bed = params.col_bed_file
        params.parent1_name = 'Col'
        log.warn "⚠️  DEPRECATED: --col_bed_file is deprecated. Use --parent1_bed and --parent1_name instead."
    }
    if (params.ler_bed_file && !params.parent2_bed) {
        params.parent2_bed = params.ler_bed_file
        params.parent2_name = 'Ler'
        log.warn "⚠️  DEPRECATED: --ler_bed_file is deprecated. Use --parent2_bed and --parent2_name instead."
    }

log.info "DEBUG: Selected profile: ${workflow.profile}"
log.info "DEBUG: Configured output dir: ${params.outdir}"
//log.info "DEBUG: Configured executor: ${workflow.config.process.executor ?: 'NOT SET'}"
//log.info "DEBUG: Configured maxForks: ${workflow.config.process.maxForks ?: 'NOT SET'}"

    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
	[],
        params.outdir,
        params.reads,
	params.sample_id,
       params.cenhapmer_db_dir,
       params.input_format   
    )

// Add these two lines immediately after the PIPELINE_INITIALISATION call
    log.info "DEBUG: PIPELINE_INITIALISATION subworkflow has finished."
    // log.info "DEBUG: Is input_tuple channel empty? ${PIPELINE_INITIALISATION.out.input_tuple.isEmpty()}" // <--- REMOVE .getVal()
    //



// --- ADD: Create combinations of ACC, CAT, CHR (now using dynamic parent names) ---
ch_kmc_combinations = Channel
    .from([params.parent1_name, params.parent2_name])
    .flatMap { acc ->
        ['CEN', 'ARMS'].collectMany { cat ->
            (1..params.num_chromosomes).collect { chr ->
                [acc, cat, chr]
            }
        }
    }
// Add debug to see what combinations are created
ch_kmc_combinations.view { "KMC combinations: $it" }
log.info "Using parent accessions: ${params.parent1_name} × ${params.parent2_name}"
log.info "Number of chromosomes: ${params.num_chromosomes}"
// --- END ADD ---



    // WORKFLOW: Run main workflow
    //
    CHARLA_NF (
        PIPELINE_INITIALISATION.out.input_tuple,
            params.cenhapmer_db_dir,
                  ch_kmc_combinations // <--- ADD this as a third input to CHARLA_NF
    )


// --- ADD THE FOLLOWING LINES ---
// Mix the two output channels from CHARLA_NF into a single channel
// This channel will contain all .kmc_pre and .kmc_suf files as they are produced
    ch_all_kmers = CHARLA_NF.out.kmc_pre.mix(CHARLA_NF.out.kmc_suf)
// --- END ADDED LINES ---

//
// SUBWORKFLOW: Run completion tasks
//
PIPELINE_COMPLETION (
    params.email,
    params.email_on_fail,
    params.plaintext_email,
    params.outdir,
    params.monochrome_logs,
    params.hook_url,
    ch_all_kmers // <--- CHANGE THIS LINE to pass the new mixed channel
)



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
