/*
    Run ARG screening tools
*/

include { ABRICATE_RUN                }  from '../../modules/nf-core/abricate/run/main'
include { AMRFINDERPLUS_UPDATE        }  from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN           }  from '../../modules/nf-core/amrfinderplus/run/main'
include { FARGENE                     }  from '../../modules/nf-core/fargene/main'
include { UNZIP                       }  from '../../modules/nf-core/unzip/main'
include { DEEPARG_DOWNLOADDATA        }  from '../../modules/nf-core/deeparg/downloaddata/main'
include { DEEPARG_PREDICT             }  from '../../modules/nf-core/deeparg/predict/main'
include { RGI_MAIN                    }  from '../../modules/nf-core/rgi/main/main'
include { HAMRONIZATION_ABRICATE      }  from '../../modules/nf-core/hamronization/abricate/main'
include { HAMRONIZATION_RGI           }  from '../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_DEEPARG       }  from '../../modules/nf-core/hamronization/deeparg/main'
include { HAMRONIZATION_AMRFINDERPLUS }  from '../../modules/nf-core/hamronization/amrfinderplus/main'
include { HAMRONIZATION_FARGENE       }  from '../../modules/nf-core/hamronization/fargene/main'
include { HAMRONIZATION_SUMMARIZE     }  from '../../modules/nf-core/hamronization/summarize/main'

workflow ARG {
    take:
    contigs // tuple val(meta), path(contigs)
    annotations // output from prokka

    main:
    ch_versions = Channel.empty()

     // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = Channel.empty()

    // AMRfinderplus run
        // Prepare channel for database
    if ( !params.arg_skip_amrfinderplus && params.arg_amrfinderplus_db ) {
        ch_amrfinderplus_db = Channel
            .fromPath( params.arg_amrfinderplus_db )
            .first()
    } else if ( !params.arg_skip_amrfinderplus && !params.arg_amrfinderplus_db ) {
        AMRFINDERPLUS_UPDATE( )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if ( !params.arg_skip_amrfinderplus ) {
        AMRFINDERPLUS_RUN ( contigs, ch_amrfinderplus_db )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

    // Reporting
        HAMRONIZATION_AMRFINDERPLUS ( AMRFINDERPLUS_RUN.out.report, 'json', AMRFINDERPLUS_RUN.out.tool_version, AMRFINDERPLUS_RUN.out.db_version )
        ch_versions = ch_versions.mix(HAMRONIZATION_AMRFINDERPLUS.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_AMRFINDERPLUS.out.json)
    }

    // fARGene run
    if ( !params.arg_skip_fargene ) {

        ch_fargene_classes = Channel.fromList( params.arg_fargene_hmmmodel.tokenize(',') )

        ch_fargene_input = contigs
                            .combine(ch_fargene_classes)
                            .map {
                                meta, contigs, hmm_class ->
                                    def meta_new = meta.clone()
                                    meta_new['hmm_class'] = hmm_class
                                [ meta_new, contigs, hmm_class ]
                            }
                            .multiMap {
                                contigs: [ it[0], it[1] ]
                                hmmclass: it[2]
                            }

        FARGENE ( ch_fargene_input.contigs, ch_fargene_input.hmmclass )
        ch_versions = ch_versions.mix(FARGENE.out.versions)

        // Reporting
        // Note: currently hardcoding versions, has to be updated with every fARGene-update
        HAMRONIZATION_FARGENE ( FARGENE.out.hmm.transpose(), 'json', '0.1', '0.1' )
        ch_versions = ch_versions.mix(HAMRONIZATION_FARGENE.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_FARGENE.out.json)
    }

    // RGI run
    if ( !params.arg_skip_rgi ) {

        RGI_MAIN ( contigs )
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

        // Reporting
        HAMRONIZATION_RGI ( RGI_MAIN.out.tsv, 'json', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version )
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_RGI.out.json)
    }

    // DeepARG prepare download
    if ( !params.arg_skip_deeparg && params.arg_deeparg_data ) {

        if ( file(params.arg_deeparg_data).getExtension() == "zip") {
            UNZIP( [ [id: "deepargdb"], params.arg_deeparg_data ] )
            ch_deeparg_db = UNZIP.out.unzipped_archive.map{meta, db -> [db]}
        } else {
            ch_deeparg_db = Channel
                .fromPath( params.arg_deeparg_data )
                .first()
        }

    } else if ( !params.arg_skip_deeparg && !params.arg_deeparg_data ) {
        DEEPARG_DOWNLOADDATA( )
        ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA.out.versions)
        ch_deeparg_db = DEEPARG_DOWNLOADDATA.out.db
    }

    // DeepARG run
    if ( !params.arg_skip_deeparg ) {

        annotations
                .map {
                    it ->
                        def meta  = it[0]
                        def anno  = it[1]
                        def model = params.arg_deeparg_model

                    [ meta, anno, model ]
                }
                .set { ch_input_for_deeparg }

        DEEPARG_PREDICT ( ch_input_for_deeparg, ch_deeparg_db )
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions)

        // Reporting
        // Note: currently hardcoding versions as unreported by DeepARG
        // Make sure to update on version bump.
        HAMRONIZATION_DEEPARG ( DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg), 'json', '1.0.2', params.arg_deeparg_data_version )
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_DEEPARG.out.json)
    }

    // ABRicate run
    if ( !params.arg_skip_abricate ) {
        ABRICATE_RUN ( contigs )
        ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)

        HAMRONIZATION_ABRICATE ( ABRICATE_RUN.out.report, 'json', '1.0.1', '2021-Mar-27' )
        ch_versions = ch_versions.mix(HAMRONIZATION_ABRICATE.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_ABRICATE.out.json)
    }

    ch_input_to_hamronization_summarize
        .map{
            it[1]
        }
        .collect()
        .set { ch_input_for_hamronization_summarize }

    HAMRONIZATION_SUMMARIZE( ch_input_for_hamronization_summarize, params.arg_hamronization_summarizeformat )
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE.out.versions)

    emit:
    versions = ch_versions
}
