/*
    Run BGC screening tools
*/

include { UNTAR as UNTAR_CSS                       } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_DETECTION                 } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_MODULES                   } from '../../modules/nf-core/modules/untar/main'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../modules/nf-core/modules/antismash/antismashlitedownloaddatabases/main'
include { ANTISMASH_ANTISMASHLITE                  } from '../../modules/nf-core/modules/antismash/antismashlite/main'
include { GECCO_RUN                                } from '../../modules/nf-core/modules/gecco/run/main'

workflow BGC {

    take:
    fna     // tuple val(meta), path(PROKKA.out.fna)
    gff     // tuple val(meta), path(PROKKA.out.gff)

    main:
    ch_versions = Channel.empty()

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

    if ( !params.bgc_skip_gecco ) {
        ch_gecco_input = fna.groupTuple()
                            .multiMap {
                                fna: [ it[0], it[1], [] ]
                            }

        GECCO_RUN ( ch_gecco_input, [] )
        ch_versions = ch_versions.mix(GECCO_RUN.out.versions)
    }

    emit:
    versions = ch_versions
}