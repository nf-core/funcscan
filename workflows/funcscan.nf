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

include { GUNZIP                  } from '../modules/nf-core/modules/gunzip/main'
include { PROKKA                  } from '../modules/nf-core/modules/prokka/main'
include { PRODIGAL                } from '../modules/nf-core/modules/prodigal/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FUNCSCAN {

    ch_versions = Channel.empty()
    ch_mqc  = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Some tools require uncompressed input
    INPUT_CHECK.out.contigs
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }
        .set { fasta_prep }

    GUNZIP ( fasta_prep.compressed )
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_prepped_input = GUNZIP.out.gunzip
                        .mix(fasta_prep.uncompressed)

    // Some tools require annotated FASTAs
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && !params.amp_skip_hmmsearch ) ) {
        if ( params.run_annotation_tool == "prodigal") {
            PRODIGAL ( ch_prepped_input, params.my_prodigal_f )
            ch_versions = ch_versions.mix(PRODIGAL.out.versions)
            ch_annotation_output = PRODIGAL.out.amino_acid_fasta
        }   else if ( params.run_annotation_tool == "prokka") {
            PROKKA ( ch_prepped_input, [], [] )
            ch_versions = ch_versions.mix(PROKKA.out.versions)
            ch_annotation_output = PROKKA.out.faa
        }
    else {
        ( ch_annotation_output = Channel.empty() )
    }
    }

    /*
        AMPs
    */
    if ( params.run_amp_screening ) {

        if ( !params.amp_skip_hmmsearch ) {
            AMP ( ch_prepped_input, ch_annotation_output )
        }   else {
            AMP ( ch_prepped_input, [] )
        }
        ch_versions = ch_versions.mix(AMP.out.versions)
        ch_mqc      = ch_mqc.mix(AMP.out.mqc)
    }

    /*
        ARGs
    */
    if ( params.run_arg_screening ) {
        if (params.arg_skip_deeparg) {
            ARG ( ch_prepped_input, [] )
        } else {
            ARG ( ch_prepped_input, ch_annotation_output )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
        ch_mqc      = ch_mqc.mix(ARG.out.mqc)
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
