/*
    Run BGC screening tools
*/

include { UNTAR as UNTAR_CSS                       } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_DETECTION                 } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_MODULES                   } from '../../modules/nf-core/untar/main'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../modules/nf-core/antismash/antismashlitedownloaddatabases/main'
include { ANTISMASH_ANTISMASHLITE                  } from '../../modules/nf-core/antismash/antismashlite/main'
include { GECCO_RUN                                } from '../../modules/nf-core/gecco/run/main'
include { HMMER_HMMSEARCH as BGC_HMMER_HMMSEARCH   } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { DEEPBGC_DOWNLOAD                         } from '../../modules/nf-core/deepbgc/download/main'
include { DEEPBGC_PIPELINE                         } from '../../modules/nf-core/deepbgc/pipeline/main'
include { COMBGC                                   } from '../../modules/local/combgc'

workflow BGC {

    take:
    fna     // tuple val(meta), path(PROKKA.out.fna)
    gff     // tuple val(meta), path(PROKKA.out.gff)
    faa     // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)

    main:
    ch_versions = Channel.empty()

    // When adding new tool that requires FAA, make sure to update conditions
    // in funcscan.nf around annotation and AMP subworkflow execution
    // to ensure annotation is executed!
    ch_faa_for_bgc_hmmsearch = faa

    // ANTISMASH
    if ( !params.bgc_skip_antismash ) {
        // Check whether user supplies database and/or antismash directory. If not, obtain them via the module antismashlite/antismashlitedownloaddatabases.
        // Important for future maintenance: For CI tests, only the "else" option below is used. Both options should be tested locally whenever the antiSMASH module gets updated.
        if ( params.bgc_antismash_databases && params.bgc_antismash_installationdirectory ) {

            ch_antismash_databases = Channel
                .fromPath( params.bgc_antismash_databases )
                .first()

            ch_antismash_directory = Channel
                .fromPath( params.bgc_antismash_installationdirectory )
                .first()

        } else {

            ch_css_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz"
            ch_detection_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz"
            ch_modules_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz"

            UNTAR_CSS ( [ [], ch_css_for_antismash ] )
            UNTAR_DETECTION ( [ [], ch_detection_for_antismash ] )
            UNTAR_MODULES ( [ [], ch_modules_for_antismash ] )

            ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR_CSS.out.untar.map{ it[1] }, UNTAR_DETECTION.out.untar.map{ it[1] }, UNTAR_MODULES.out.untar.map{ it[1] } )
            ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.versions)
            ch_antismash_databases = ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database

            ch_antismash_directory = ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir

        }

        ch_antismash_input = fna.mix(gff)
                                .groupTuple()
                                .multiMap {
                                    fna: [ it[0], it[1][0] ]
                                    gff: it[1][1]
                                }

        ANTISMASH_ANTISMASHLITE ( ch_antismash_input.fna, ch_antismash_databases, ch_antismash_directory, ch_antismash_input.gff )
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITE.out.versions)
    }

    // DEEPBGC
    if ( !params.bgc_skip_deepbgc ){
        if ( params.bgc_deepbgc_database ) {

            ch_deepbgc_database = Channel
                .fromPath( params.bgc_deepbgc_database )
                .first()
        } else {
            DEEPBGC_DOWNLOAD()
            ch_deepbgc_database = DEEPBGC_DOWNLOAD.out.db
        }

    ch_deepbgc_input = fna

    DEEPBGC_PIPELINE ( ch_deepbgc_input, ch_deepbgc_database)
    ch_versions = ch_versions.mix(DEEPBGC_PIPELINE.out.versions)
    }

    // GECCO
    if ( !params.bgc_skip_gecco ) {
        ch_gecco_input = fna.groupTuple()
                            .multiMap {
                                fna: [ it[0], it[1], [] ]
                            }

        GECCO_RUN ( ch_gecco_input, [] )
        ch_versions = ch_versions.mix(GECCO_RUN.out.versions)
    }

    // HMMSEARCH
    if ( !params.bgc_skip_hmmsearch ) {
        if ( params.bgc_hmmsearch_models ) { ch_bgc_hmm_models = Channel.fromPath( params.bgc_hmmsearch_models, checkIfExists: true ) } else { exit 1, '[nf-core/funcscan] error: hmm model files not found for --bgc_hmmsearch_models! Please check input.' }

        ch_bgc_hmm_models_meta = ch_bgc_hmm_models
            .map {
                file ->
                    def meta  = [:]
                    meta['id'] = file.extension == 'gz' ? file.name - '.hmm.gz' :  file.name - '.hmm'

                [ meta, file ]
            }

        ch_in_for_bgc_hmmsearch = ch_faa_for_bgc_hmmsearch.combine(ch_bgc_hmm_models_meta)
            .map {
                meta_faa, faa, meta_hmm, hmm ->
                    def meta_new = [:]
                    meta_new['id']     = meta_faa['id']
                    meta_new['hmm_id'] = meta_hmm['id']
                // TODO make optional outputs params?
                [ meta_new, hmm, faa, params.bgc_hmmsearch_savealignments, params.bgc_hmmsearch_savetargets, params.bgc_hmmsearch_savedomains ]
            }

        BGC_HMMER_HMMSEARCH ( ch_in_for_bgc_hmmsearch )
        ch_versions = ch_versions.mix(BGC_HMMER_HMMSEARCH.out.versions)
    }

    // COMBGC
    if ( !params.bgc_skip_antismash ) { ch_antismash = Channel.fromPath("${params.outdir}/bgc/antismash/", type: 'dir') } else { ch_antismash = [] }
    if ( !params.bgc_skip_deepbgc )   { ch_deepbgc   = Channel.fromPath("${params.outdir}/bgc/deepbgc/", type: 'dir') } else { ch_deepbgc = [] }
    if ( !params.bgc_skip_gecco )     { ch_gecco     = Channel.fromPath("${params.outdir}/bgc/gecco/", type: 'dir') } else { ch_gecco = [] }
    ch_outdir = Channel.fromPath("${params.outdir}/reports/combgc/", type: 'dir')

    COMBGC( ch_antismash, ch_deepbgc, ch_gecco, ch_outdir )

    emit:
    versions = ch_versions
}
