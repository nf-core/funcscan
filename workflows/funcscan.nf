/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFuncscan.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { AMP } from '../subworkflows/local/amp'
include { ARG } from '../subworkflows/local/arg'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { GUNZIP as GUNZIP_FASTA      } from '../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_FAA        } from '../modules/nf-core/modules/gunzip/main'
include { PROKKA                      } from '../modules/nf-core/modules/prokka/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FUNCSCAN {

    ch_versions = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    /* TODO Separate out FASTA & FAA using multiMAP
        -> separte GUNZIP for FAA
        -> downstream re-merge with .combine
        -> conditional logic which to run PROKKA/which NOT
        -> ?multiMap input into AMP workflows? Or change expected input tuple?

    */

    // Some tools require uncompressed input. Here separate out compressed and
    // uncompressed files accordingly for FASTA/FAA then mix them all back
    // together again to a single channel for downstream.
    ch_prep_for_gunzip = INPUT_CHECK.out.contigs
        .multiMap {
            meta, fasta, faa ->
                fasta: [ meta, fasta ]
                faa: [ meta, faa ]
        }

    ch_fasta_for_gunzip = ch_prep_for_gunzip.fasta
        .branch {
            meta, file ->
                compressed: file.extension == 'gz'
                uncompressed: true
        }

    ch_faa_for_gunzip = ch_prep_for_gunzip.faa
        .branch {
            meta, file ->
                compressed: file != '' && file.extension == 'gz'
                uncompressed: true
        }

    GUNZIP_FASTA ( ch_fasta_for_gunzip.compressed )
    GUNZIP_FAA ( ch_faa_for_gunzip.compressed )
    ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    ch_versions = ch_versions.mix(GUNZIP_FAA.out.versions)

    ch_prepped_input_fasta = GUNZIP_FASTA.out.gunzip.mix(ch_fasta_for_gunzip.uncompressed)
    ch_prepped_input_faa = GUNZIP_FAA.out.gunzip.mix(ch_faa_for_gunzip.uncompressed)
    ch_unzipped_input = ch_prepped_input_fasta.join(ch_prepped_input_faa)

    // Some tools require annotated FASTAs, but Prokka is slow - so only run
    // when tools that use them activated
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && !params.amp_skip_hmmsearch ) ) {

        // Only run if we haven't already got annotated files already
        ch_input_to_prokka = ch_unzipped_input
            .branch{
                meta, fasta, faa ->
                    run: faa == ''
                    skip: true
            }

        // Extract just the fasta to match PROKKA expected input
        PROKKA ( ch_input_to_prokka.run.map {meta, fasta, faa -> [meta, fasta]}, [], [] )

        // Join back newly created FAA files with the original contigs for
        // downstream
        ch_output_from_prokka = ch_input_to_prokka.run
            .map {
                meta, fasta, faa ->
                    [ meta, fasta ]
            }
            .join( PROKKA.out.faa )
        ch_prepped_input = ch_input_to_prokka.skip.mix( ch_output_from_prokka )
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt)
    } else {
        ch_prepped_input = ch_unzipped_input
    }

    /*
        AMPs
    */
    if ( params.run_amp_screening ) {

        if ( !params.amp_skip_hmmsearch ) {
            AMP ( ch_prepped_input )
        } else {
            AMP ( ch_prepped_input, [] )
        }


        ch_versions = ch_versions.mix(AMP.out.versions)
        ch_multiqc_files      = ch_multiqc_files.mix(AMP.out.mqc)
    }

    /*
        ARGs
    */
    // TODO WHY DO WE NEED PROKKA.OUT.FNA? HOW TO SEND DOWNSTREAM? WHAT HAPPENS
    // IF PRE-SUPPLIED FAA, SO PROKKA NOT RUN? APPARENTLY ONLY NEEDS
    if ( params.run_arg_screening ) {
        ARG ( ch_prepped_input, PROKKA.out.fna )
        ch_versions = ch_versions.mix(ARG.out.versions)
        ch_multiqc_files      = ch_multiqc_files.mix(ARG.out.mqc)
    }

    /*
        BGCs
    */
    // TODO antismash

    // Cleaning up versions
    CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.unique().collectFile(name: 'collated_versions.yml') )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFuncscan.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
