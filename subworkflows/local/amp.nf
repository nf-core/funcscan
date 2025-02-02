/*
    Run AMP screening tools
*/

include { MACREL_CONTIGS                                              } from '../../modules/nf-core/macrel/contigs/main'
include { HMMER_HMMSEARCH as AMP_HMMER_HMMSEARCH                      } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { AMPLIFY_PREDICT                                             } from '../../modules/nf-core/amplify/predict/main'
include { AMPIR                                                       } from '../../modules/nf-core/ampir/main'
include { AMP_DATABASE_DOWNLOAD                                       } from '../../modules/local/amp_database_download'
include { AMPCOMBI2_PARSETABLES                                       } from '../../modules/nf-core/ampcombi2/parsetables'
include { AMPCOMBI2_COMPLETE                                          } from '../../modules/nf-core/ampcombi2/complete'
include { AMPCOMBI2_CLUSTER                                           } from '../../modules/nf-core/ampcombi2/cluster'
include { GUNZIP as GUNZIP_MACREL_PRED ; GUNZIP as GUNZIP_MACREL_ORFS } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as AMP_GUNZIP_HMMER_HMMSEARCH                        } from '../../modules/nf-core/gunzip/main'
include { TABIX_BGZIP as AMP_TABIX_BGZIP                              } from '../../modules/nf-core/tabix/bgzip/main'
include { MERGE_TAXONOMY_AMPCOMBI                                     } from '../../modules/local/merge_taxonomy_ampcombi'

workflow AMP {
    take:
    fastas          // tuple val(meta), path(contigs)
    faas            // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)
    tsvs            // tuple val(meta), path(MMSEQS_CREATETSV.out.tsv)
    gbks            // tuple val(meta), path(ANNOTATION_ANNOTATION_TOOL.out.gbk)
    tsvs_interpro   // tuple val(meta), path(INTERPROSCAN.out.tsv)'

    main:
    ch_versions                    = Channel.empty()
    ch_ampresults_for_ampcombi     = Channel.empty()
    ch_macrel_faa                  = Channel.empty()
    ch_ampcombi_summaries          = Channel.empty()
    ch_ampcombi_complete           = null

    // When adding new tool that requires FAA, make sure to update conditions
    // in funcscan.nf around annotation and AMP subworkflow execution
    // to ensure annotation is executed!
    ch_faa_for_amplify             = faas
    ch_faa_for_amp_hmmsearch       = faas
    ch_faa_for_ampir               = faas
    ch_faa_for_ampcombi            = faas
    ch_gbk_for_ampcombi            = gbks
    ch_interpro_for_ampcombi       = tsvs_interpro

    // AMPLIFY
    if ( !params.amp_skip_amplify ) {
        AMPLIFY_PREDICT ( ch_faa_for_amplify, [] )
        ch_versions                = ch_versions.mix( AMPLIFY_PREDICT.out.versions )
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix( AMPLIFY_PREDICT.out.tsv )
    }

    // MACREL
    if ( !params.amp_skip_macrel ) {
        MACREL_CONTIGS ( fastas )
        ch_versions                = ch_versions.mix( MACREL_CONTIGS.out.versions )
        GUNZIP_MACREL_PRED ( MACREL_CONTIGS.out.amp_prediction )
        GUNZIP_MACREL_ORFS ( MACREL_CONTIGS.out.all_orfs )
        ch_versions                = ch_versions.mix( GUNZIP_MACREL_PRED.out.versions )
        ch_versions                = ch_versions.mix( GUNZIP_MACREL_ORFS.out.versions )
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix( GUNZIP_MACREL_PRED.out.gunzip )
        ch_macrel_faa              = ch_macrel_faa.mix( GUNZIP_MACREL_ORFS.out.gunzip )
        ch_faa_for_ampcombi        = ch_faa_for_ampcombi.mix( ch_macrel_faa )
    }

    // AMPIR
    if ( !params.amp_skip_ampir ) {
        AMPIR ( ch_faa_for_ampir, params.amp_ampir_model, params.amp_ampir_minlength, 0.0 )
        ch_versions                = ch_versions.mix( AMPIR.out.versions )
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix( AMPIR.out.amps_tsv )
    }

    // HMMSEARCH
    if ( params.amp_run_hmmsearch ) {
        if ( params.amp_hmmsearch_models ) { ch_amp_hmm_models = Channel.fromPath( params.amp_hmmsearch_models, checkIfExists: true ) } else { error('[nf-core/funcscan] error: HMM model files not found for --amp_hmmsearch_models! Please check input.') }

        ch_amp_hmm_models_meta = ch_amp_hmm_models
            .map {
                file ->
                    def meta   = [:]
                    meta['id'] = file.extension == 'gz' ? file.name - '.hmm.gz' :  file.name - '.hmm'
                [ meta, file ]
            }

        ch_in_for_amp_hmmsearch = ch_faa_for_amp_hmmsearch
                                    .combine( ch_amp_hmm_models_meta )
                                    .map {
                                        meta_faa, faa, meta_hmm, hmm ->
                                            def meta_new = [:]
                                            meta_new['id']     = meta_faa['id']
                                            meta_new['hmm_id'] = meta_hmm['id']
                                        [ meta_new, hmm, faa, params.amp_hmmsearch_savealignments, params.amp_hmmsearch_savetargets, params.amp_hmmsearch_savedomains ]
                                    }

        AMP_HMMER_HMMSEARCH ( ch_in_for_amp_hmmsearch )
        ch_versions = ch_versions.mix( AMP_HMMER_HMMSEARCH.out.versions )
        AMP_GUNZIP_HMMER_HMMSEARCH ( AMP_HMMER_HMMSEARCH.out.output )
        ch_versions = ch_versions.mix( AMP_GUNZIP_HMMER_HMMSEARCH.out.versions )
        ch_AMP_GUNZIP_HMMER_HMMSEARCH = AMP_GUNZIP_HMMER_HMMSEARCH.out.gunzip
            .map { meta, file ->
                [ [id: meta.id], file ]
            }
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix( ch_AMP_GUNZIP_HMMER_HMMSEARCH )
    }

    // AMPCOMBI2
    ch_input_for_ampcombi = ch_ampresults_for_ampcombi
        .groupTuple()
        .join( ch_faa_for_ampcombi )
        .join( ch_gbk_for_ampcombi )
        .join( ch_interpro_for_ampcombi )
        .multiMap{
            input: [ it[0], it[1] ]
            faa: it[2]
            gbk: it[3]
            interpro: it [4]
        }

    // AMPCOMBI2::PARSETABLES
    if ( params.amp_ampcombi_db != null ) {
        AMPCOMBI2_PARSETABLES ( ch_input_for_ampcombi.input,  ch_input_for_ampcombi.faa,  ch_input_for_ampcombi.gbk, params.amp_ampcombi_db_id, params.amp_ampcombi_db, ch_input_for_ampcombi.interpro )
        } else {
            AMP_DATABASE_DOWNLOAD( params.amp_ampcombi_db_id )
            ch_versions = ch_versions.mix( AMP_DATABASE_DOWNLOAD.out.versions )
            ch_ampcombi_input_db = AMP_DATABASE_DOWNLOAD.out.db
            AMPCOMBI2_PARSETABLES ( ch_input_for_ampcombi.input, ch_input_for_ampcombi.faa, ch_input_for_ampcombi.gbk, params.amp_ampcombi_db_id, ch_ampcombi_input_db, ch_input_for_ampcombi.interpro )
        }
    ch_versions = ch_versions.mix( AMPCOMBI2_PARSETABLES.out.versions )

    ch_ampcombi_summaries = AMPCOMBI2_PARSETABLES.out.tsv.map{ it[1] }.collect()

    // AMPCOMBI2::COMPLETE
    ch_summary_count = ch_ampcombi_summaries.map { it.size() }.sum()

    if ( ch_summary_count == 0 || ch_summary_count == 1 )  {
        log.warn("[nf-core/funcscan] AMPCOMBI2: ${ch_summary_count} file(s) passed. Skipping AMPCOMBI2_COMPLETE, AMPCOMBI2_CLUSTER, and TAXONOMY MERGING steps.")
    } else {
        AMPCOMBI2_COMPLETE(ch_ampcombi_summaries)
        ch_versions = ch_versions.mix( AMPCOMBI2_COMPLETE.out.versions )
        ch_ampcombi_complete = AMPCOMBI2_COMPLETE.out.tsv
                                .filter { file -> file.countLines() > 1 }
    }

    // AMPCOMBI2::CLUSTER
    if ( ch_ampcombi_complete != null )  {
        AMPCOMBI2_CLUSTER ( ch_ampcombi_complete )
        ch_versions = ch_versions.mix( AMPCOMBI2_CLUSTER.out.versions )
    } else {
        log.warn("[nf-core/funcscan] No AMP hits were found in the samples and so no clustering will be applied.")
    }

    // MERGE_TAXONOMY
    if ( params.run_taxa_classification && ch_ampcombi_complete == null ) {
        log.warn("[nf-core/funcscan] No AMP hits were found in the samples, therefore no Taxonomy will be merged ")
    } else if ( params.run_taxa_classification && ch_ampcombi_complete != null ) {
        ch_mmseqs_taxonomy_list = tsvs.map{ it[1] }.collect()

        MERGE_TAXONOMY_AMPCOMBI( AMPCOMBI2_CLUSTER.out.cluster_tsv, ch_mmseqs_taxonomy_list )
        ch_versions = ch_versions.mix( MERGE_TAXONOMY_AMPCOMBI.out.versions )

        ch_tabix_input = Channel.of( [ 'id':'ampcombi_complete_summary_taxonomy' ] )
            .combine( MERGE_TAXONOMY_AMPCOMBI.out.tsv )

        AMP_TABIX_BGZIP( ch_tabix_input )
        ch_versions = ch_versions.mix( AMP_TABIX_BGZIP.out.versions )
    }

    emit:
    versions = ch_versions
}
