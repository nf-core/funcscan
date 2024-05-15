/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_funcscan_pipeline'
include { paramsSummaryMap; validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
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
include { TAXA_CLASS } from '../subworkflows/local/taxa_class'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_INPUT_PREP } from '../modules/nf-core/gunzip/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_LONG  } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_SHORT } from '../modules/nf-core/seqkit/seq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FUNCSCAN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_input = Channel.fromSamplesheet("input")

    ///////////////////////
    // INPUT PREPARATION //
    ///////////////////////

    // Some tools require uncompressed input
    ch_input_prep = ch_input
                        .map{meta, fasta, faa, gbk -> [meta, [fasta, faa, gbk]]}
                        .transpose()
                        .branch {
                            compressed: it[1].toString().endsWith('.gz')
                            uncompressed: it[1]
                        }

    GUNZIP_INPUT_PREP ( ch_input_prep.compressed )
    ch_versions = ch_versions.mix( GUNZIP_INPUT_PREP.out.versions )

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_intermediate_input = GUNZIP_INPUT_PREP.out.gunzip
                            .mix( ch_input_prep.uncompressed )
                            .groupTuple()
                            .map{
                                meta, files ->
                                    def fasta_found   = files.find{it.toString().tokenize('.').last().matches('fasta|fas|fna|fa')}
                                    def faa_found     = files.find{it.toString().endsWith('.faa')}
                                    def gbk_found     = files.find{it.toString().tokenize('.').last().matches('gbk|gbff')}
                                    def fasta         = fasta_found   != null ? fasta_found   : []
                                    def faa           = faa_found     != null ? faa_found     : []
                                    def gbk           = gbk_found     != null ? gbk_found     : []

                                    [meta, fasta, faa, gbk]
                            }
                            .branch {
                                meta, fasta, faa, gbk ->
                                    preannotated: gbk != []
                                    fastas: true
                            }

    // Split each FASTA into long and short contigs to
    // speed up BGC workflow with BGC-compatible contig lengths only
    ch_intermediate_fasta_for_split = ch_intermediate_input.fastas.map{ meta, fasta, faa, gbk -> [ meta, fasta ] }
    SEQKIT_SEQ_LONG ( ch_intermediate_fasta_for_split )
    SEQKIT_SEQ_SHORT ( ch_intermediate_fasta_for_split )
    ch_versions = ch_versions.mix( SEQKIT_SEQ_LONG.out.versions )
    ch_versions = ch_versions.mix( SEQKIT_SEQ_SHORT.out.versions )

    ch_intermediate_input_long = SEQKIT_SEQ_LONG.out.fastx
                                .map{ meta, file -> [ meta + [id: meta.id + '_long', length: "long" ], file ] }
                                .filter{
                                    meta, fasta ->
                                        !fasta.isEmpty()
                                }

    ch_intermediate_input_short = SEQKIT_SEQ_SHORT.out.fastx
                                .map{ meta, file -> [ meta + [id: meta.id + '_short', length: "short" ], file ] }
                                .filter{
                                    meta, fasta ->
                                        !fasta.isEmpty()
                                }

    // Now they are split, can annotated together for efficiency
    ch_input_for_annotation = ch_intermediate_input_long.mix( ch_intermediate_input_short )

    /*
        ANNOTATION
    */

    // Some tools require annotated FASTAs
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {

        ANNOTATION( ch_input_for_annotation )
        ch_versions = ch_versions.mix( ANNOTATION.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( ANNOTATION.out.multiqc_files )

        ch_new_annotation = ch_input_for_annotation
                                .join( ANNOTATION.out.faa )
                                .join( ANNOTATION.out.gbk )

    } else {
        ch_new_annotation = Channel.empty()
    }

    // Mix back the preannotated samples with the newly annotated ones,
    // but also have dedicated channel for subworkflows that should only use
    // for long contigs
    ch_prepped_input = ch_intermediate_input.preannotated
                        .mix( ch_new_annotation )
                        .multiMap {
                            meta, fasta, faa, gbk ->
                                fastas: [meta, fasta]
                                faas: [meta, faa]
                                gbks: [meta, gbk]
                        }

    ch_prepped_input_long = ch_new_annotation
                                .filter{
                                    meta, fasta, faa, gbk ->
                                        meta.length == "long"
                                }
                                .mix(ch_intermediate_input.preannotated)
                                .map {
                                    meta, fasta, faa, gbk ->
                                        if ( params.run_bgc_screening && meta.length == null ) { log.warn("[nf-core/funcscan] Pre-annotated input will not be filtered to long contigs for BGC screening! Expect long-run times and/or possible crashes if includes very short contigs") }
                                    [meta, fasta, faa, gbk]
                                }
                                .multiMap {
                                    meta, fasta, faa, gbk ->
                                        fastas: [meta, fasta]
                                        faas: [meta, faa]
                                        gbks: [meta, gbk]
                                }

    /*
        TAXONOMIC CLASSIFICATION
    */

    // The final subworkflow reports need taxonomic classification.
    // This can be either on NT or AA level depending on annotation.
    if ( params.run_taxa_classification ) {

            if ( params.run_bgc_screening && !params.run_amp_screening && !params.run_arg_screening ) {
                ch_input_for_taxonomy = ch_prepped_input_long.fastas
            } else {
                ch_input_for_taxonomy = ch_prepped_input.fastas
            }

            TAXA_CLASS ( ch_input_for_taxonomy )
            ch_versions     = ch_versions.mix( TAXA_CLASS.out.versions )
            ch_taxonomy_tsv = TAXA_CLASS.out.sample_taxonomy

    } else {
            ch_mmseqs_db              = Channel.empty()
            ch_taxonomy_querydb       = Channel.empty()
            ch_taxonomy_querydb_taxdb = Channel.empty()
            ch_taxonomy_tsv           = Channel.empty()
    }

    ///////////////
    // SCREENING //
    ///////////////

    /*
        AMPs
    */
    if ( params.run_amp_screening && !params.run_taxa_classification ) {
        AMP (
            ch_prepped_input.fastas,
            ch_prepped_input.faas
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()

                },
            ch_taxonomy_tsv
        )
        ch_versions = ch_versions.mix(AMP.out.versions)
    } else if ( params.run_amp_screening && params.run_taxa_classification ) {
        AMP (
            ch_prepped_input.fastas,
            ch_prepped_input.faas
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                    },
            ch_taxonomy_tsv
                .filter {
                        meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                    }
        )
        ch_versions = ch_versions.mix( AMP.out.versions )
    }

    /*
        ARGs
    */
    if ( params.run_arg_screening && !params.run_taxa_classification ) {
        if ( params.arg_skip_deeparg ) {
            ARG (
                ch_prepped_input.fastas,
                [],
                ch_taxonomy_tsv
                )
        } else {
            ARG (
                ch_prepped_input.fastas,
                ch_prepped_input.faas
                    .filter {
                        meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. ARG screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    },
                ch_taxonomy_tsv
            )
        }
        ch_versions = ch_versions.mix( ARG.out.versions )
    } else if ( params.run_arg_screening && params.run_taxa_classification ) {
        if ( params.arg_skip_deeparg ) {
            ARG (
                ch_prepped_input.fastas,
                [],
                ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                    }
                )
        } else {
            ARG (
                ch_prepped_input.fastas,
                ch_prepped_input.faas
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. ARG screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    },
                ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                }
            )
        }
        ch_versions = ch_versions.mix( ARG.out.versions )
    }

    /*
        BGCs
    */
    if ( params.run_bgc_screening && !params.run_taxa_classification ) {
        BGC (
            ch_prepped_input_long.fastas,
            ch_prepped_input_long.faas
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_prepped_input_long.gbks
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty GBK file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_taxonomy_tsv
        )
        ch_versions = ch_versions.mix( BGC.out.versions )
    } else if ( params.run_bgc_screening && params.run_taxa_classification ) {
        BGC (
            ch_prepped_input.fastas,
            ch_prepped_input.faas
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced produced an empty FAA file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_prepped_input.gbks
                .filter {
                    meta, file ->
                        if ( file != [] && file.isEmpty() ) log.warn("[nf-core/funcscan] Annotation of following sample produced an empty GBK file. BGC screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                },
            ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("[nf-core/funcscan] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                }
        )
        ch_versions = ch_versions.mix( BGC.out.versions )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML( ch_versions )
        .collectFile( storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //

    ch_multiqc_config                     = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    summary_params                        = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json" )
    ch_workflow_summary                   = Channel.value( paramsSummaryMultiqc( summary_params ) )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value( methodsDescriptionText( ch_multiqc_custom_methods_description ))
    ch_multiqc_files                      = ch_multiqc_files.mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
    ch_multiqc_files                      = ch_multiqc_files.mix( ch_collated_versions )
    ch_multiqc_files                      = ch_multiqc_files.mix( ch_methods_description.collectFile(name: 'methods_description_mqc.yaml') )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
