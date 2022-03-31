/*
    Run AMP screening tools
*/

include { MACREL_CONTIGS          } from '../../modules/nf-core/modules/macrel/contigs/main'

workflow AMP {
    take:
    contigs // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()
    ch_mqc      = Channel.empty()

    // TODO AMPEP(?)
    // TODO ampir
    if ( !params.amp_skip_macrel ) {
        MACREL_CONTIGS ( contigs )
        ch_versions = ch_versions.mix(MACREL_CONTIGS.out.versions)
    }

    emit:
    versions = ch_versions
    mqc = ch_mqc

}
