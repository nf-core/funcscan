/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFuncscan.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Validate fARGene inputs
// Split input into array, find the union with our valid classes, extract only
// invalid classes, and if they exist, exit. Note `tokenize` used here as this
// works for `interesect` and other groovy functions, but require `split` for
// `Channel.of` creation. See `arg.nf` for latter.
def fargene_classes = params.arg_fargene_hmmmodel
def fargene_valid_classes = [ "class_a", "class_b_1_2", "class_b_3",
                            "class_c", "class_d_1", "class_d_2",
                            "qnr", "tet_efflux", "tet_rpg", "tet_enzyme"
                            ]
def fargene_user_classes = fargene_classes.tokenize(',')
def fargene_classes_valid = fargene_user_classes.intersect( fargene_valid_classes )
def fargene_classes_missing = fargene_user_classes - fargene_classes_valid

if ( fargene_classes_missing.size() > 0 ) exit 1, "[nf-core/funcscan] ERROR: invalid class present in --arg_fargene_hmmodel. Please check input. Invalid class: ${fargene_classes_missing.join(', ')}"


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
include { BGC } from '../subworkflows/local/bgc'

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
    ch_funcscan_logo = Channel.fromPath("$projectDir/docs/images/nf-core-funcscan_logo_flat_light.png")

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

    /*
        ANNOTATION
    */

    // Some tools require annotated FASTAs
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening ) ) {

        if ( params.run_annotation_tool == "prodigal") {
            PRODIGAL ( ch_prepped_input, "gff" )
            ch_versions              = ch_versions.mix(PRODIGAL.out.versions)
            ch_annotation_faa        = PRODIGAL.out.amino_acid_fasta
            ch_annotation_fna        = PRODIGAL.out.nucleotide_fasta
            ch_annotation_gff        = PRODIGAL.out.gene_annotations
        }   else if ( params.run_annotation_tool == "prokka") {
            PROKKA ( ch_prepped_input, [], [] )
            ch_versions              = ch_versions.mix(PROKKA.out.versions)
            ch_annotation_faa        = PROKKA.out.faa
            ch_annotation_fna        = PROKKA.out.fna
            ch_annotation_gff        = PROKKA.out.gff
        }

    } else {

        ch_annotation_faa        = Channel.empty()
        ch_annotation_fna        = Channel.empty()
        ch_annotation_gff        = Channel.empty()

    }

    /*
        SCREENING
    */

    /*
        AMPs
    */
    if ( params.run_amp_screening ) {
        if ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) {
            AMP ( ch_prepped_input, ch_annotation_faa )
        }   else {
            AMP ( ch_prepped_input, [] )
        }
        ch_versions = ch_versions.mix(AMP.out.versions)
    }

    /*
        ARGs
    */
    if ( params.run_arg_screening ) {
        if (params.arg_skip_deeparg) {
            ARG ( ch_prepped_input, [] )
        } else {
            ARG ( ch_prepped_input, ch_annotation_faa )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
    }

    /*
        BGCs
    */
    if ( params.run_bgc_screening ) {
        BGC ( fasta_prep.compressed, ch_annotation_gff )
        ch_version = ch_versions.mix(BGC.out.versions)
    }

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

    ch_multiqc_files = ch_multiqc_files.mix(ch_funcscan_logo)

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
