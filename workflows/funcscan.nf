/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

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
include { ANNOTATION } from '../subworkflows/local/annotation'
include { AMP        } from '../subworkflows/local/amp'
include { ARG        } from '../subworkflows/local/arg'
include { BGC        } from '../subworkflows/local/bgc'

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
include { GUNZIP as GUNZIP_INPUT_PREP       } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FNA     } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FAA     } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_GFF     } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FNA    } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FAA    } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_GFF    } from '../modules/nf-core/gunzip/main'
include { BIOAWK                            } from '../modules/nf-core/bioawk/main'
include { PROKKA                            } from '../modules/nf-core/prokka/main'
include { PRODIGAL as PRODIGAL_GFF          } from '../modules/nf-core/prodigal/main'
include { PRODIGAL as PRODIGAL_GBK          } from '../modules/nf-core/prodigal/main'
include { PYRODIGAL                         } from '../modules/nf-core/pyrodigal/main'
include { BAKTA_BAKTADBDOWNLOAD             } from '../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA                       } from '../modules/nf-core/bakta/bakta/main'

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

    ch_input = Channel.fromSamplesheet("input")

    ///////////////////////
    // INPUT PREPARATION //
    ///////////////////////

    // Some tools require uncompressed input
    ch_input_prep = ch_input
                        .map{meta, fasta, faa, feature -> [meta, [fasta, faa, feature]]}
                        .transpose()
                        .branch {
                            compressed: it[1].toString().endsWith('.gz')
                            uncompressed: it[1]
                        }

    GUNZIP_INPUT_PREP ( ch_input_prep.compressed )
    ch_versions = ch_versions.mix(GUNZIP_INPUT_PREP.out.versions)

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_intermediate_input = GUNZIP_INPUT_PREP.out.gunzip
                            .mix(ch_input_prep.uncompressed)
                            .groupTuple()
                            .map{
                                meta, files ->
                                    def fasta_found   = files.find{it.toString().tokenize('.').last().matches('fasta|fas|fna|fa')}
                                    def faa_found     = files.find{it.toString().endsWith('.faa')}
                                    def feature_found = files.find{it.toString().tokenize('.').last().matches('gff|gbk')}

                                    // https://github.com/antismash/antismash/issues/364
                                    if ( params.run_bgc_screening && !params.bgc_skip_antismash && feature_found != null ) {
                                        log.warn("[nf-core/funcscan] antiSMASH screening requested and pre-annotated files given.")
                                        log.warn("Be aware that Prokka generated GFF or GBK files will likely fail with antiSMASH!")
                                        log.warn("See usage docs. File: " + feature_found.name) }

                                    def fasta   = fasta_found   != null ? fasta_found : []
                                    def faa     = faa_found     != null ? faa_found : []
                                    def feature = feature_found != null ? feature_found : []

                                    [meta, fasta, faa, feature]
                            }
                            .multiMap {
                                meta, fasta, faa, feature ->
                                    fastas: [ meta, fasta ]
                                    annotations : [ meta, faa, feature ]
                            }

    // Add to meta the length of longest contig for downstream filtering
    ch_intermediate_input.fastas
    ch_intermediate_input.annotations

    BIOAWK ( ch_intermediate_input.fastas )
    ch_versions = ch_versions.mix(BIOAWK.out.versions)

    ch_intermediate_input = ch_intermediate_input.fastas
                                .join(BIOAWK.out.longest)
                                .join(ch_intermediate_input.annotations)
                                .map{
                                    meta, fasta, length, faa, feature ->
                                        def meta_new = [:]
                                        meta_new['longest_contig'] = Integer.parseInt(length)
                                    [ meta + meta_new, fasta, faa, feature ]
                                }

    ////////////////
    // ANNOTATION //
    ////////////////

    // Separate pre-annotated FASTAs from those that need annotation
    ch_input_for_annotation = ch_intermediate_input
                                .branch {
                                    meta, fasta, protein, feature ->
                                        preannotated: protein != []
                                        unannotated: true
                                }

    // Some tools require annotated FASTAs
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {

        ch_unannotated_for_annotation = ch_input_for_annotation.unannotated
                                            .map{
                                                meta, fasta, protein, feature ->
                                                [meta, fasta]
                                            }

        ANNOTATION( ch_unannotated_for_annotation )
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)

        // Only Bakta and Prokka make GBK, else give empty entry to satisfy downstream cardinality
        if ( ['bakta', 'prokka'].contains(params.annotation_tool) ) {
            ch_new_annotation = ch_unannotated_for_annotation
                                    .join(ANNOTATION.out.faa)
                                    .join(ANNOTATION.out.gff)
                                    .join(ANNOTATION.out.gbk)
        } else {
            ch_new_annotation = ch_unannotated_for_annotation
                        .join(ANNOTATION.out.faa)
                        .join(ANNOTATION.out.gff)
                        .map {
                            meta, fasta, faa, gff ->
                                [meta, fasta, faa, gff, []]
                        }
        }

    } else {
        ch_new_annotation = Channel.empty()
    }

    ch_prepped_input = ch_input_for_annotation.preannotated
                        .map{
                            meta, fasta, protein, feature ->
                                def gff = feature.extension == 'gff' ? feature : []
                                def gbk = feature.extension == 'gbk' ? feature : []
                            [meta, fasta, protein, gff, gbk]
                        }
                        .mix(ch_new_annotation)
                        .multiMap {
                            meta, fasta, protein, gff, gbk ->
                            fastas: [meta, fasta]
                            faas: [meta, protein]
                            gffs: [meta, gff]
                            gbks: [meta, gbk]
                        }

    ///////////////
    // SCREENING //
    ///////////////

    /*
        AMPs
    */

    if ( params.run_amp_screening ) {
        AMP (
            ch_prepped_input.fastas,
            ch_prepped_input.faas
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
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
            ARG ( ch_prepped_input.fastas, [] )
        } else {
            ARG (
                ch_prepped_input.fastas,
                ch_prepped_input.faas
                    .filter {
                        meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. ARG screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    }
            )
        }
        ch_versions = ch_versions.mix(ARG.out.versions)
    }

    // /*
    //     BGCs
    // */
    if ( params.run_bgc_screening ) {
        BGC (
            ch_prepped_input.fastas,
            ch_prepped_input.faas
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_prepped_input.gffs
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty GFF file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_prepped_input.gbks
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty GBK file. AMP screening tools requiring this file will not be executed: ${meta.id}")
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


    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {

        if ( params.annotation_tool == 'prokka' ) {
            ch_multiqc_files = ch_multiqc_files.mix(ANNOTATION.out.multiqc_files.map{it[1]})
        }

    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {
        if( ['prokka','bakta'].contains(params.annotation_tool) ){
            ch_multiqc_files = ch_multiqc_files.mix( ANNOTATION.out.multiqc_files.collect{it[1]}.ifEmpty([]))
        }
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
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
