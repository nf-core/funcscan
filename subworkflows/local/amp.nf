/*
    Run AMP screening tools
*/

include { MACREL_CONTIGS          } from '../../modules/nf-core/modules/macrel/contigs/main'
include { HMMER_HMMSEARCH         } from '../../modules/nf-core/modules/hmmer/hmmsearch/main'
include { AMPLIFY_PREDICT         } from '../../modules/nf-core/modules/amplify/predict/main'
include { AMPIR                   } from '../../modules/nf-core/modules/ampir/main'

workflow AMP {
    take:
    contigs // tuple val(meta), path(contigs)
    faa     // tuple val(meta), path(PROKKA.out.faa)

    main:
    ch_versions = Channel.empty()
    ch_mqc      = Channel.empty()

    ch_faa_for_amplify = faa
    ch_faa_for_hmmsearch = faa
    ch_faa_for_ampir = faa

    // AMPLIFY
    if ( !params.amp_skip_amplify ) {
        AMPLIFY_PREDICT ( ch_faa_for_amplify, [] )
        ch_versions = ch_versions.mix(AMPLIFY_PREDICT.out.versions)
    }

    // MACREL
    if ( !params.amp_skip_macrel ) {
        MACREL_CONTIGS ( contigs )
        ch_versions = ch_versions.mix(MACREL_CONTIGS.out.versions)
    }

    // HMMSEARCH
    if ( !params.amp_skip_hmmsearch ) {
        if (params.amp_hmmsearch_models) { ch_amp_hmm_models = Channel.fromPath( params.amp_hmmsearch_models, checkIfExists: true ) } else { exit 1, '[nf-core/funscan] error: hmm model files not found for --amp_hmmsearch_models! Please check input.' }

        ch_amp_hmm_models_meta = ch_amp_hmm_models
            .map {
                file ->
                    def meta  = [:]
                    meta['id'] = file.extension == 'gz' ? file.name - '.hmm.gz' :  file.name - '.hmm'

                [ meta, file ]
            }

        ch_in_for_hmmsearch = ch_faa_for_hmmsearch.combine(ch_amp_hmm_models_meta)
            .map {
                meta_faa, faa, meta_hmm, hmm ->
                    def meta_new = [:]
                    meta_new['id'] = meta_faa['id'] + '_'  + meta_hmm['id']
                // TODO make optional outputs params?
                [ meta_new, hmm, faa, params.amp_hmmsearch_savealignments, params.amp_hmmsearch_savetargets, params.amp_hmmsearch_savedomains ]
            }

        HMMER_HMMSEARCH ( ch_in_for_hmmsearch )
    }

    // AMPIR
    if ( !params.amp_skip_ampir ) {
        ch_faa_for_ampir.dump(tag: "faa_for_ampir")
        println(params.amp_ampir_minlength)
        println(params.amp_ampir_minprobability)
        println(params.amp_ampir_model)
        //AMPIR ( ch_faa_for_ampir, params.amp_ampir_model, params.amp_ampir_minlength, params.amp_ampir_minprobability )
    }

    emit:
    versions = ch_versions
    mqc = ch_mqc

}
