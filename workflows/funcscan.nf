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
include { AMP        }  from '../subworkflows/local/amp'
include { ARG        }  from '../subworkflows/local/arg'
include { BGC        }  from '../subworkflows/local/bgc'
include { TAXA_CLASS } from '../subworkflows/local/taxa_class'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { GUNZIP as GUNZIP_FASTA_PREP    } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FNA  } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_FAA  } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRODIGAL_GFF  } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FNA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_FAA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PYRODIGAL_GFF } from '../modules/nf-core/gunzip/main'
include { PROKKA                         } from '../modules/nf-core/prokka/main'
include { PRODIGAL as PRODIGAL_GFF       } from '../modules/nf-core/prodigal/main'
include { PRODIGAL as PRODIGAL_GBK       } from '../modules/nf-core/prodigal/main'
include { PYRODIGAL as PYRODIGAL_GBK     } from '../modules/nf-core/pyrodigal/main'
include { PYRODIGAL as PYRODIGAL_GFF     } from '../modules/nf-core/pyrodigal/main'
include { BAKTA_BAKTADBDOWNLOAD          } from '../modules/nf-core/bakta/baktadbdownload/main'
include { BAKTA_BAKTA                    } from '../modules/nf-core/bakta/bakta/main'
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

    // Some tools require uncompressed input
    fasta_prep = ch_input
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }

    GUNZIP_FASTA_PREP ( fasta_prep.compressed )
    ch_versions = ch_versions.mix( GUNZIP_FASTA_PREP.out.versions )

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_unzipped_fastas = GUNZIP_FASTA_PREP.out.gunzip
                        .mix( fasta_prep.uncompressed )

    // Split each FASTA into long and short contigs to
    // speed up BGC workflow with BGC-compatible contig lengths only
    SEQKIT_SEQ_LONG ( ch_unzipped_fastas )
    SEQKIT_SEQ_SHORT ( ch_unzipped_fastas )
    ch_versions = ch_versions.mix( SEQKIT_SEQ_LONG.out.versions )
    ch_versions = ch_versions.mix( SEQKIT_SEQ_SHORT.out.versions )

    ch_prepped_input_long = SEQKIT_SEQ_LONG.out.fastx
                                .map{ meta, file -> [ meta + [id: meta.id + '_long', length: "long" ], file ] }
                                .filter{
                                    meta, fasta ->
                                        !fasta.isEmpty()
                                }

    ch_prepped_input_short = SEQKIT_SEQ_SHORT.out.fastx
                                .map{ meta, file -> [ meta + [id: meta.id + '_short', length: "short" ], file ]}
                                .filter{
                                    meta, fasta ->
                                        !fasta.isEmpty()
                                }

    ch_prepped_input = ch_prepped_input_long.mix( ch_prepped_input_short )

    /*
        TAXONOMIC CLASSIFICATION
    */

    // The final subworkflow reports need taxonomic classification.
    // This can be either on NT or AA level depending on annotation.
    // TODO: Only NT at the moment. AA tax. classification will be added only when its PR is merged.
    if ( params.run_taxa_classification ) {
            TAXA_CLASS ( ch_prepped_input )
            ch_versions     = ch_versions.mix( TAXA_CLASS.out.versions )
            ch_taxonomy_tsv = TAXA_CLASS.out.sample_taxonomy

    } else {

            ch_mmseqs_db              = Channel.empty()
            ch_taxonomy_querydb       = Channel.empty()
            ch_taxonomy_querydb_taxdb = Channel.empty()
            ch_taxonomy_tsv           = Channel.empty()
    }

    /*
        ANNOTATION
    */

    // Some tools require annotated FASTAs
    // For prodigal: run twice, once for gff and once for gbk generation, (for parity with PROKKA which produces both)
    if ( ( params.run_arg_screening && !params.arg_skip_deeparg ) || ( params.run_amp_screening && ( !params.amp_skip_hmmsearch || !params.amp_skip_amplify || !params.amp_skip_ampir ) ) || ( params.run_bgc_screening && ( !params.bgc_skip_hmmsearch || !params.bgc_skip_antismash ) ) ) {

        if ( params.annotation_tool == "prodigal" ) {
            PRODIGAL_GFF ( ch_prepped_input, "gff" )
            GUNZIP_PRODIGAL_FAA ( PRODIGAL_GFF.out.amino_acid_fasta )
            GUNZIP_PRODIGAL_FNA ( PRODIGAL_GFF.out.nucleotide_fasta )
            GUNZIP_PRODIGAL_GFF ( PRODIGAL_GFF.out.gene_annotations )
            ch_versions              = ch_versions.mix( PRODIGAL_GFF.out.versions )
            ch_annotation_faa        = GUNZIP_PRODIGAL_FAA.out.gunzip
            ch_annotation_fna        = GUNZIP_PRODIGAL_FNA.out.gunzip
            ch_annotation_gff        = GUNZIP_PRODIGAL_GFF.out.gunzip
            ch_annotation_gbk        = Channel.empty() // Prodigal GBK and GFF output are mutually exclusive

            if ( params.save_annotations == true ) {
                PRODIGAL_GBK ( ch_prepped_input, "gbk" )
                ch_versions          = ch_versions.mix( PRODIGAL_GBK.out.versions )
                ch_annotation_gbk    = PRODIGAL_GBK.out.gene_annotations // Prodigal GBK output stays zipped because it is currently not used by any downstream subworkflow.
            }
        } else if ( params.annotation_tool == "pyrodigal" ) {
            PYRODIGAL_GFF ( ch_prepped_input, "gff" )
            GUNZIP_PYRODIGAL_FAA ( PYRODIGAL_GFF.out.faa )
            GUNZIP_PYRODIGAL_FNA ( PYRODIGAL_GFF.out.fna )
            GUNZIP_PYRODIGAL_GFF ( PYRODIGAL_GFF.out.annotations )
            ch_versions              = ch_versions.mix( PYRODIGAL_GFF.out.versions )
            ch_annotation_faa        = GUNZIP_PYRODIGAL_FAA.out.gunzip
            ch_annotation_fna        = GUNZIP_PYRODIGAL_FNA.out.gunzip
            ch_annotation_gff        = GUNZIP_PYRODIGAL_GFF.out.gunzip
            ch_annotation_gbk        = Channel.empty() // Pyrodigal GBK and GFF output are mutually exclusive

            if ( params.save_annotations == true ) {
                PYRODIGAL_GBK ( ch_prepped_input, "gbk" )
                ch_versions          = ch_versions.mix( PYRODIGAL_GBK.out.versions )
                ch_annotation_gbk    = PYRODIGAL_GBK.out.annotations // Pyrodigal GBK output stays zipped because it is currently not used by any downstream subworkflow.
            }
        }  else if ( params.annotation_tool == "prokka" ) {
            PROKKA ( ch_prepped_input, [], [] )
            ch_versions              = ch_versions.mix( PROKKA.out.versions )
            ch_annotation_faa        = PROKKA.out.faa
            ch_annotation_fna        = PROKKA.out.fna
            ch_annotation_gff        = PROKKA.out.gff
            ch_annotation_gbk        = PROKKA.out.gbk
        }   else if ( params.annotation_tool == "bakta" ) {

            // BAKTA prepare download
            if ( params.annotation_bakta_db_localpath ) {
                ch_bakta_db          = Channel
                    .fromPath( params.annotation_bakta_db_localpath )
                    .first()
            } else {
                BAKTA_BAKTADBDOWNLOAD ( )
                ch_versions          = ch_versions.mix( BAKTA_BAKTADBDOWNLOAD.out.versions )
                ch_bakta_db          = ( BAKTA_BAKTADBDOWNLOAD.out.db )
            }

            BAKTA_BAKTA ( ch_prepped_input, ch_bakta_db, [], [] )
            ch_versions              = ch_versions.mix( BAKTA_BAKTA.out.versions )
            ch_annotation_faa        = BAKTA_BAKTA.out.faa
            ch_annotation_fna        = BAKTA_BAKTA.out.fna
            ch_annotation_gff        = BAKTA_BAKTA.out.gff
            ch_annotation_gbk        = BAKTA_BAKTA.out.gbff
        }

    } else {

        ch_annotation_faa            = Channel.empty()
        ch_annotation_fna            = Channel.empty()
        ch_annotation_gff            = Channel.empty()
        ch_annotation_gbk            = Channel.empty()

    }

    /*
        SCREENING
    */

    /*
        AMPs
    */
    if ( params.run_amp_screening && !params.run_taxa_classification ) {
        AMP (
            ch_prepped_input,
            ch_annotation_faa
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()

                },
            ch_taxonomy_tsv
        )
        ch_versions = ch_versions.mix(AMP.out.versions)
    } else if ( params.run_amp_screening && params.run_taxa_classification ) {
        AMP (
            ch_prepped_input,
            ch_annotation_faa
                .filter {
                    meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                        !file.isEmpty()
                    },
            ch_taxonomy_tsv
                .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
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
                ch_prepped_input,
                [],
                ch_taxonomy_tsv
                )
        } else {
            ARG (
                ch_prepped_input,
                ch_annotation_faa
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    },
                ch_taxonomy_tsv
            )
        }
        ch_versions = ch_versions.mix( ARG.out.versions )
    } else if ( params.run_arg_screening && params.run_taxa_classification ) {
        if ( params.arg_skip_deeparg ) {
            ARG (
                ch_prepped_input,
                [],
                ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                    }
                )
        } else {
            ARG (
                ch_prepped_input,
                ch_annotation_faa
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Annotation of following sample produced produced an empty FAA file. AMP screening tools requiring this file will not be executed: ${meta.id}")
                            !file.isEmpty()
                    },
                ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
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
            ch_prepped_input_long,
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
                },
            ch_taxonomy_tsv
        )
        ch_versions = ch_versions.mix( BGC.out.versions )
    } else if ( params.run_bgc_screening && params.run_taxa_classification ) {
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
                },
            ch_taxonomy_tsv
                    .filter {
                        meta, file ->
                        if ( file.isEmpty() ) log.warn("Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                        !file.isEmpty()
                }
        )
        ch_versions = ch_versions.mix( BGC.out.versions )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("${workflow.projectDir}/docs/images/nf-core-funcscan_logo_light.png", checkIfExists: true)

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    if( params.annotation_tool=='prokka' ) {
        ch_multiqc_files                  = ch_multiqc_files.mix( PROKKA.out.txt.collect{it[1]}.ifEmpty([]) )
    }

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
