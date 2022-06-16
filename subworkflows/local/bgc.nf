/*
    Run AMP screening tools
*/
include { UNTAR as UNTAR_CSS                       } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_DETECTION                 } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_MODULES                   } from '../../modules/nf-core/modules/untar/main'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../modules/nf-core/modules/antismash/antismashlitedownloaddatabases/main'
include { ANTISMASH_ANTISMASHLITE                  } from '../../modules/nf-core/modules/antismash/antismashlite/main'

workflow BGC {

    take:
    fna     // tuple val(meta), path(PROKKA.out.fna)
    gff     // tuple val(meta), path(PROKKA.out.gff)

    main:
    ch_versions = Channel.empty()
    ch_mqc      = Channel.empty()
    // Check whether user supplies database and/or antismash directory. If not, use those from nf-core.
    // params.bgc_antismash_databases = null
    // params.bgc_antismash_directory = null
    ch_antismash_databases = Channel
        .fromPath( params.bgc_antismash_databases )
        .first()
    ch_antismash_directory = Channel
        .fromPath( params.bgc_antismash_directory )
        .first()

    if ( !ch_antismash_databases ) {

        ch_css_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz"
        ch_detection_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz"
        ch_modules_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz"

        UNTAR_CSS ( [ [], ch_css_for_antismash ] )
        UNTAR_DETECTION ( [ [], ch_detection_for_antismash ] )
        UNTAR_MODULES ( [ [], ch_modules_for_antismash ] )

        ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR_CSS.out.untar.map{ it[1] }, UNTAR_DETECTION.out.untar.map{ it[1] }, UNTAR_MODULES.out.untar.map{ it[1] } )
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.versions)
        ch_antismash_databases = ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database

        if ( !ch_antismash_directory) {
            ch_antismash_directory = ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir
        }
    }

    ANTISMASH_ANTISMASHLITE ( fna, ch_antismash_databases, ch_antismash_directory, gff.map{ it[1] } )
        ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITE.out.versions)

    emit:
    versions = ch_versions
    mqc = ch_mqc

}
