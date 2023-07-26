/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowFuncscan.initialise(params, log)

// Check input path parameters to see if they exist
/*def checkPathParamList = [ params.input, params.multiqc_config, params.annotation_bakta_db_localpath,
                            params.amp_hmmsearch_models, params.amp_ampcombi_db,
                            params.arg_amrfinderplus_db, params.arg_deeparg_data,
                            params.bgc_antismash_databases, params.bgc_antismash_installationdirectory,
                            params.bgc_deepbgc_database, params.bgc_hmmsearch_models ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { error("Input samplesheet not specified!") }
*/

// Validate fARGene inputs
// Split input into array, find the union with our valid classes, extract only
// invalid classes, and if they exist, exit. Note `tokenize` used here as this
// works for `interesect` and other groovy functions, but require `split` for
// `Channel.of` creation. See `arg.nf` for latter.
/*
def fargene_classes = params.arg_fargene_hmmmodel
def fargene_valid_classes = [ "class_a", "class_b_1_2", "class_b_3",
                            "class_c", "class_d_1", "class_d_2",
                            "qnr", "tet_efflux", "tet_rpg", "tet_enzyme"
                            ]
def fargene_user_classes = fargene_classes.tokenize(',')
def fargene_classes_valid = fargene_user_classes.intersect( fargene_valid_classes )
def fargene_classes_missing = fargene_user_classes - fargene_classes_valid

if ( fargene_classes_missing.size() > 0 ) error("[nf-core/funcscan] ERROR: invalid class present in --arg_fargene_hmmodel. Please check input. Invalid class: ${fargene_classes_missing.join(', ')}")
*/

// Validate antiSMASH inputs
// 1. Make sure that either both or none of the antiSMASH directories are supplied
if ( ( params.run_bgc_screening && !params.bgc_antismash_databases && params.bgc_antismash_installationdirectory && !params.bgc_skip_antismash) || ( params.run_bgc_screening && params.bgc_antismash_databases && !params.bgc_antismash_installationdirectory && !params.bgc_skip_antismash ) ) error("[nf-core/funcscan] ERROR: You supplied either the antiSMASH database or its installation directory, but not both. Please either supply both directories or none (letting the pipeline download them instead).")

// 2. If both are supplied: Exit if we have a name collision error
else if ( params.run_bgc_screening && params.bgc_antismash_databases && params.bgc_antismash_installationdirectory && !params.bgc_skip_antismash ) {
    antismash_database_dir = new File(params.bgc_antismash_databases)
    antismash_install_dir = new File(params.bgc_antismash_installationdirectory)
    if ( antismash_database_dir.name == antismash_install_dir.name ) error("[nf-core/funcscan] ERROR: Your supplied antiSMASH database and installation directories have identical names: \"" + antismash_install_dir.name + "\".\nPlease make sure to name them differently, for example:\n - Database directory:      "+ antismash_database_dir.parent + "/antismash_db\n - Installation directory:  " + antismash_install_dir.parent + "/antismash_dir")
}

// 3. Give warning if not using container system assuming conda

if ( params.run_bgc_screening && ( !params.bgc_antismash_databases || !params.bgc_antismash_installationdirectory ) && !params.bgc_skip_antismash && ( session.config.conda && session.config.conda.enabled ) ) { log.warn "[nf-core/funcscan] Running antiSMASH download database module, and detected conda has been enabled. Assuming using conda for pipeline run, check config if this is not expected!" }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ANNOTATION                    } from '../subworkflows/local/annotation'
include { AMP                           } from '../subworkflows/local/amp'
include { ARG                           } from '../subworkflows/local/arg'
include { BGC                           } from '../subworkflows/local/bgc'
include { GUNZIP as GUNZIP_FASTA_PREP   } from '../modules/nf-core/gunzip/main'
include { BIOAWK                        } from '../modules/nf-core/bioawk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FUNCSCAN {

    ch_versions = Channel.empty()
    ch_multiqc_logo = Channel.fromPath("$projectDir/docs/images/nf-core-funcscan_logo_flat_light.png")

    /*
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
*/
    ch_input = Channel.fromSamplesheet("input")

    // Some tools require uncompressed input
    ch_fasta_prep = ch_input
        .map {
            meta, fasta, protein, feature ->
            [meta, fasta]
        }
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }

    GUNZIP_FASTA_PREP ( ch_fasta_prep.compressed )
    ch_versions = ch_versions.mix(GUNZIP_FASTA_PREP.out.versions)

    ch_preannotated_files = ch_input
                                .filter {
                                    meta, fasta, protein, feature ->
                                    (protein || feature)
                                }
                                .map {
                                    meta, fasta, protein, feature ->
                                    [meta, protein, feature]
                                }
                                .dump(tag: "ch_preannotated_files")

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_prepped_fastas = GUNZIP_FASTA_PREP.out.gunzip
                        .mix(ch_fasta_prep.uncompressed)

    // Add to meta the length of longest contig for downstream filtering
    BIOAWK ( ch_prepped_fastas )
    ch_versions = ch_versions.mix(BIOAWK.out.versions)

    ch_prepped_input = ch_prepped_fastas
                        .join( BIOAWK.out.longest )
                        .map{
                            meta, fasta, length ->
                                def meta_new = meta.clone()
                                meta['longest_contig'] = Integer.parseInt(length)
                            [ meta, fasta ]
                        }
                        .dump(tag: "ch_prepped_input")

    /*
        ANNOTATION
    */
    // TODO NOT JOINING HERE BECAUSE OF BIOAWK MISSING IN INFO CH_PREANNOTATED_FILES
    // Separate out already annotated from unannotated
    ch_input_for_annotation = ch_prepped_input // or original ch_prepped_fastas, but fastas missing BIOAWK
                                .join(ch_preannotated_files)
                                .dump(tag: "ch_input_for_annotation_prebranch")
                                .branch {
                                    meta, fasta, protein, feature ->
                                        annotated_protein: protein != []
                                        annotated_feature: feature != []
                                        unannotated: true
                                }

    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {

        ANNOTATION ( ch_input_for_annotation.unannotated.dump(tag: "ch_input_for_annotation.unannotated_premap").map{meta, fasta, protein, feature -> [meta, fasta]} )
        ch_new_annotation_faa   = ANNOTATION.out.faa
        ch_new_annotation_fna   = ANNOTATION.out.fna
        ch_new_annotation_gff   = ANNOTATION.out.gff
        ch_new_annotation_gbk   = ANNOTATION.out.gbk
        ch_versions             = ANNOTATION.out.versions

    } else {

        ch_new_annotation_faa   = Channel.empty()
        ch_new_annotation_fna   = Channel.empty()
        ch_new_annotation_gff   = Channel.empty()
        ch_new_annotation_gbk   = Channel.empty()

    }

    // Join back the pre-annotated FASTAs with newly annotated FASTAs
    ch_annotation_proteins = ch_input_for_annotation.annotated_feature.map{ meta, fasta, protein, feature -> [meta, protein] }
    ch_annotation_faa      = ch_new_annotation_faa.mix(ch_annotation_proteins)

    ch_annotation_features = ch_input_for_annotation.annotated_feature.map{ meta, fasta, protein, feature -> [meta, feature] }
    ch_annotation_gff      = ch_annotation_features.filter { meta, feature -> feature.toString().endsWith('.gff') }.mix(ch_new_annotation_gff)
    ch_annotation_gbk      = ch_annotation_features.filter { meta, feature -> feature.toString().endsWith('.gbk') }.mix(ch_new_annotation_gbk)

    /*
        SCREENING
    */

    /*
        AMPs
    */
    if ( params.run_amp_screening ) {
        AMP (
            ch_prepped_input,
            ch_annotation_faa
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                }
        )
        ch_versions = ch_versions.mix(AMP.out.versions)
    }

    /*
        ARGs
    */
    if ( params.run_arg_screening ) {
        if (params.arg_skip_deeparg) {
            ARG ( ch_prepped_input, [] )
        } else {
            ARG (
                ch_prepped_input,
                ch_annotation_faa
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    }
            )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
    }

    /*
        BGCs
    */
    if ( params.run_bgc_screening ) {
        BGC (
            ch_prepped_input,
            ch_annotation_gff
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty GFF file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_annotation_faa
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_annotation_gbk
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty GBK file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                }
        )
        ch_versions = ch_versions.mix(BGC.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFuncscan.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowFuncscan.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {
        ch_multiqc_files = ch_multiqc_files.mix( ANNOTATION.out.mqc.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
