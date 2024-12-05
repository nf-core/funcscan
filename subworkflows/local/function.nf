/*
    RUN FUNCTIONAL CLASSIFICATION
*/

include { INTERPROSCAN  } from '../../modules/nf-core/interproscan/main'

workflow FUNCTION {
    take:
    faas // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)

    main:
    ch_versions                 = Channel.empty()
    ch_interproscan_tsv         = Channel.empty()

    ch_faa_for_interproscan     = faas

    INTERPROSCAN( ch_faa_for_interproscan, [] )

    ch_interproscan_tsv = ch_interproscan_tsv.mix(INTERPROSCAN.out.tsv)

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    emit:
    versions         = ch_versions
    interproscan_tsv = ch_interproscan_tsv // channel: [ val(meta), tsv ]
}
