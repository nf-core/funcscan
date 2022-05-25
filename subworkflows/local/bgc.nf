/*
    Run AMP screening tools
*/
include { UNTAR as UNTAR1 } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR2 } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR3 } from '../../modules/nf-core/modules/untar/main'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../modules/nf-core/modules/antismash/antismashlitedownloaddatabases/main'
include { ANTISMASH_ANTISMASHLITE } from '../../modules/nf-core/modules/antismash/antismashlite/main'

workflow BGC {
    take:
    fna     // tuple val(meta), path(PROKKA.out.fna)
    gff     // tuple val(meta), path(PROKKA.out.gff)

    main:
    ch_versions = Channel.empty()
    ch_mqc      = Channel.empty()

    ch_css_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz"
    ch_detection_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz"
    ch_modules_for_antismash = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz"
    ch_fna_for_antismash = fna
    ch_gff_for_antismash = gff

    UNTAR1 ( [ [], ch_css_for_antismash ] )
    UNTAR2 ( [ [], ch_detection_for_antismash ] )
    UNTAR3 ( [ [], ch_modules_for_antismash ] )

    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] }, UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.versions)

    ANTISMASH_ANTISMASHLITE ( ch_fna_for_antismash, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir, ch_gff_for_antismash.map{ it[1] } )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITE.out.versions)

    emit:
    versions = ch_versions
    mqc = ch_mqc

}
