/*
    Run ARG screening tools
*/

include { FARGENE                 } from '../../modules/nf-core/modules/fargene/main'
include { DEEPARG_DOWNLOADDATA    } from '../../modules/nf-core/modules/deeparg/downloaddata/main'
include { DEEPARG_PREDICT         } from '../../modules/nf-core/modules/deeparg/predict/main'

include { HAMRONIZATION_DEEPARG   } from '../../modules/nf-core/modules/hamronization/deeparg/main'
include { HAMRONIZATION_SUMMARIZE } from '../../modules/nf-core/modules/hamronization/summarize/main'

workflow ARG {
    take:
    contigs // file: /path/to/samplesheet.csv
    annotations // output from prokka

    main:
    ch_versions = Channel.empty()
    ch_mqc      = Channel.empty()

     // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = Channel.empty()

    // fARGene run
    if ( !params.skip_arg_fargene ) {
        FARGENE ( contigs, params.fargene_hmm_model )
        ch_versions = ch_versions.mix(FARGENE.out.versions)
    }

    // DeepARG prepare download
    if ( !params.skip_arg_deeparg && params.deeparg_data ) {
        Channel
            .fromPath( params.deeparg_data )
            .set { ch_deeparg_db }
    } else if ( !params.skip_arg_deeparg && !params.deeparg_data ) {
        DEEPARG_DOWNLOADDATA( )
        DEEPARG_DOWNLOADDATA.out.db.set { ch_deeparg_db }
    }

    // DeepARG run

    if ( !params.skip_arg_deeparg ) {

        annotations
                .map {
                    it ->
                        def meta  = it[0]
                        def anno  = it[1]
                        def model = params.deeparg_model

                    [ meta, anno, model ]
                }
                .set { ch_input_for_deeparg }

        DEEPARG_PREDICT ( ch_input_for_deeparg, ch_deeparg_db )
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions)

    // Reporting
    // Note:currently hardcoding versions
    // how to automate in the future - but DEEPARG won't change as abandonware?
        HAMRONIZATION_DEEPARG ( DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg).dump(tag: "in_hamr_deep"), 'json', '1.0.2', '2'  )
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_DEEPARG.out.json)
    }

    ch_input_to_hamronization_summarize
        .dump(tag: "map_in")
        .map{
            it[1]
        }
        .collect()
        .dump(tag: "map_out")
        .set { ch_input_for_hamronization_summarize }

    HAMRONIZATION_SUMMARIZE( ch_input_for_hamronization_summarize, params.hamronization_summarize_format )

    emit:
    versions = ch_versions
    mqc = ch_mqc

}
