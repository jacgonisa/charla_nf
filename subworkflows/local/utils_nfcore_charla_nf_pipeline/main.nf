//
// Subworkflow with functionality specific to the jacgonisa/charla_nf pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    args
    outdir            //  string: The output directory where the results will be saved
    reads               //  string: Path to input reads (e.g., fastq file)
    sample_id           //  string: Sample identifier
    cenhapmer_db_dir
    input_format 
 
   main:
   // log.error "!!! FORCED CRASH DEBUG: This message should definitely appear in the log file. !!!"
   // throw new Exception("!!! INTENTIONAL CRASH FOR DEBUGGING PURPOSES. !!!")
    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
//    UTILS_NEXTFLOW_PIPELINE (
//        version,
//        true,
//        outdir,
//        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
//    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
//    UTILS_NFSCHEMA_PLUGIN (
//        workflow,
//        validate_params,
//        null
//    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    



    
    log.info "DEBUG: Inside PIPELINE_INITIALISATION subworkflow"
    log.info "DEBUG: 'reads' parameter value: ${reads}"
    log.info "DEBUG: Does file exist at this path? ${file(reads).exists()}"
//    log.info "DEBUG: Full path attempted: ${file(reads).getAbsolutePath()}"

    //
    // Create channel from input file provided through params.input
    //

    emit:
    input_tuple = Channel.of(tuple(file(reads), sample_id, params.kmer_size))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//

// Custom validation for pipeline parameters
//
def validateInputParameters() {

    // Validate params.reads
    if (!params.reads) {
        error "No input reads file specified! Please provide with --reads."
    } else {
        def readsFile = file(params.reads)
        if (!readsFile.exists()) {
            error "Provided reads file '${params.reads}' does not exist! Please check the path."
        }
        if (!readsFile.canRead()) {
            error "Provided reads file '${params.reads}' is not readable! Please check file permissions."
        }
    }

    // Validate params.sample_id
    if (!params.sample_id) {
        error "No sample identifier specified! Please provide with --sample_id."
    } else if (params.sample_id.trim() == '') {
        error "Sample identifier '--sample_id' cannot be empty."
    }

    // Validate params.kmer_size
    if (!params.kmer_size) {
        error "No k-mer size specified! Please provide with --kmer_size."
    } else {
        try {
            def kmerSizeInt = params.kmer_size as int
            if (kmerSizeInt <= 0) {
                error "K-mer size '--kmer_size' must be a positive integer."
            }
            // Optional: K-mer size is often odd for many tools (e.g. Jellyfish)
            if (kmerSizeInt % 2 == 0) {
                log.warn "K-mer size '--kmer_size' is an even number. Some k-mer tools prefer odd k-mer sizes. Consider using an odd number."
            }
        } catch (MissingMethodException | IllegalArgumentException e) {
            error "Invalid k-mer size '--kmer_size': Must be an integer."
        }
    }
}



//
// Validate channels from input samplesheet
//
//def validateInputSamplesheet(input) {
//    def (metas, fastqs) = input[1..2]

//    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
//    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
//    if (!endedness_ok) {
//        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
//    }

//    return [ metas[0], fastqs ]
//}
//
