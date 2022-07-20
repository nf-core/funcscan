/*
    Run ARG screening tools
*/

include { FARGENE                                  } from '../../modules/nf-core/modules/fargene/main'
include { DEEPARG_DOWNLOADDATA                     } from '../../modules/nf-core/modules/deeparg/downloaddata/main'
include { DEEPARG_PREDICT                          } from '../../modules/nf-core/modules/deeparg/predict/main'
include { RGI_MAIN                                 } from '../../modules/nf-core/modules/rgi/main/main'
include { HAMRONIZATION_RGI                        } from '../../modules/nf-core/modules/hamronization/rgi/main'
include { HAMRONIZATION_DEEPARG                    } from '../../modules/nf-core/modules/hamronization/deeparg/main'
include { HAMRONIZATION_SUMMARIZE                  } from '../../modules/nf-core/modules/hamronization/summarize/main'

workflow ARG {
    take:
    contigs // tuple val(meta), path(contigs)
    annotations // output from prokka

    main:
    ch_versions = Channel.empty()

     // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = Channel.empty()

    // fARGene run
    if ( !params.arg_skip_fargene ) {

        ch_fargene_classes = Channel.fromList( params.arg_fargene_hmmmodel.tokenize(',') )

        ch_fargene_input = contigs
                            .dump(tag: "fargene_contigs_raw")
                            .combine(ch_fargene_classes)
                            .dump(tag: "fargene_contigs_hmmclass")
                            .map {
                                meta, contigs, hmm_class ->
                                    def meta_new = meta.clone()
                                    meta_new['hmm_class'] = hmm_class
                                [ meta_new, contigs, hmm_class ]
                            }
                            .dump(tag: "fargene_updated_meta")
                            .multiMap {
                                contigs: [ it[0], it[1] ]
                                hmmclass: it[2]
                            }

        FARGENE ( ch_fargene_input.contigs, ch_fargene_input.hmmclass )
        ch_versions = ch_versions.mix(FARGENE.out.versions)

    }

    // RGI run
    if ( !params.arg_skip_rgi ) {

        RGI_MAIN ( contigs )
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

    // Reporting
    // Note: currently hardcoding versions, has to be updated with every RGI-Container-update
        HAMRONIZATION_RGI ( RGI_MAIN.out.tsv, 'json', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version )
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_RGI.out.json)

    }

    // DeepARG prepare download
    if ( !params.arg_skip_deeparg && params.arg_deeparg_data ) {
        ch_deeparg_db = Channel
            .fromPath( params.arg_deeparg_data )
            .first()
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
    // Note:currently hardcoding versions
    // how to automate in the future - but DEEPARG won't change as abandonware?
        HAMRONIZATION_DEEPARG ( DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg), 'json', '1.0.2', '2'  )
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_DEEPARG.out.json)
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
