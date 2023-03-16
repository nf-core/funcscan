/*
    Run AMP screening tools
*/

include { MACREL_CONTIGS                                            } from '../../modules/nf-core/macrel/contigs/main'
include { HMMER_HMMSEARCH as AMP_HMMER_HMMSEARCH                    } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { AMPLIFY_PREDICT                                           } from '../../modules/nf-core/amplify/predict/main'
include { AMPIR                                                     } from '../../modules/nf-core/ampir/main'
include { DRAMP_DOWNLOAD                                            } from '../../modules/local/dramp_download'
include { AMPCOMBI                                                  } from '../../modules/nf-core/ampcombi/main'
include { GUNZIP as GUNZIP_MACREL ; GUNZIP as GUNZIP_HMMER          } from '../../modules/nf-core/gunzip/main'

workflow AMP {
    take:
    contigs // tuple val(meta), path(contigs)
    faa     // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)

    main:
    ch_versions                    = Channel.empty()
    ch_ampresults_for_ampcombi     = Channel.empty()
    ch_ampcombi_summaries          = Channel.empty()

    // When adding new tool that requires FAA, make sure to update conditions
    // in funcscan.nf around annotation and AMP subworkflow execution
    // to ensure annotation is executed!
    ch_faa_for_amplify          = faa
    ch_faa_for_amp_hmmsearch    = faa
    ch_faa_for_ampir            = faa
    ch_faa_for_ampcombi         = faa

    // AMPLIFY
    if ( !params.amp_skip_amplify ) {
        AMPLIFY_PREDICT ( ch_faa_for_amplify, [] )
        ch_versions = ch_versions.mix(AMPLIFY_PREDICT.out.versions)
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix(AMPLIFY_PREDICT.out.tsv)
    }

    // MACREL
    if ( !params.amp_skip_macrel ) {
        MACREL_CONTIGS ( contigs )
        ch_versions = ch_versions.mix(MACREL_CONTIGS.out.versions)
        GUNZIP_MACREL ( MACREL_CONTIGS.out.amp_prediction )
        ch_versions = ch_versions.mix(GUNZIP_MACREL.out.versions)
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix(GUNZIP_MACREL.out.gunzip)
    }

    // AMPIR
    if ( !params.amp_skip_ampir ) {
        AMPIR ( ch_faa_for_ampir, params.amp_ampir_model, params.amp_ampir_minlength, 0.0 )
        ch_versions = ch_versions.mix(AMPIR.out.versions)
        ch_ampresults_for_ampcombi = ch_ampresults_for_ampcombi.mix(AMPIR.out.amps_tsv)
    }

    // HMMSEARCH
    if ( !params.amp_skip_hmmsearch ) {
        if ( params.amp_hmmsearch_models ) { ch_amp_hmm_models = Channel.fromPath( params.amp_hmmsearch_models, checkIfExists: true ) } else { exit 1, '[nf-core/funcscan] error: hmm model files not found for --amp_hmmsearch_models! Please check input.' }

        ch_amp_hmm_models_meta = ch_amp_hmm_models
            .map {
                file ->
                    def meta  = [:]
                    meta['id'] = file.extension == 'gz' ? file.name - '.hmm.gz' :  file.name - '.hmm'

                [ meta, file ]
            }

        ch_in_for_amp_hmmsearch = ch_faa_for_amp_hmmsearch.combine(ch_amp_hmm_models_meta)
            .map {
                meta_faa, faa, meta_hmm, hmm ->
                    def meta_new = [:]
                    meta_new['id']     = meta_faa['id']
                    meta_new['hmm_id'] = meta_hmm['id']
                [ meta_new, hmm, faa, params.amp_hmmsearch_savealignments, params.amp_hmmsearch_savetargets, params.amp_hmmsearch_savedomains ]
            }

        AMP_HMMER_HMMSEARCH ( ch_in_for_amp_hmmsearch )
        ch_versions = ch_versions.mix(AMP_HMMER_HMMSEARCH.out.versions)
    }

    //AMPCOMBI
    ch_input_for_ampcombi = ch_ampresults_for_ampcombi
        .groupTuple()
        .join( ch_faa_for_ampcombi )
        .multiMap{
            input: [ it[0], it[1] ]
            faa: it[2]
        }
    // Checks if `--amp_database` is a user supplied path and if the path does not exist it goes to default, which downloads the DRAMP database once.
    if ( params.amp_ampcombi_db ) {
        ch_ampcombi_input_db = Channel
                                    .fromPath( params.amp_ampcombi_db, checkIfExists: true ) }
    else {
        DRAMP_DOWNLOAD()
        ch_ampcombi_input_db = DRAMP_DOWNLOAD.out.db
    }

    AMPCOMBI( ch_input_for_ampcombi.input, ch_input_for_ampcombi.faa, ch_ampcombi_input_db )
    ch_ampcombi_summaries = ch_ampcombi_summaries.mix(AMPCOMBI.out.csv)

    //AMPCOMBI concatenation
    ch_ampcombi_summaries_out = ch_ampcombi_summaries
        .multiMap{
                input: [ it[0] ]
                summary: it[1]
            }
    ch_ampcombi_summaries_out.summary.collectFile(name: 'ampcombi_complete_summary.csv', storeDir: "${params.outdir}/reports/ampcombi", keepHeader:true)

    emit:
    versions = ch_versions

}
