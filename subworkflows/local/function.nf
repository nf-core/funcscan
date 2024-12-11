/*
    RUN FUNCTIONAL CLASSIFICATION
*/

include { INTERPROSCAN }          from '../../modules/nf-core/interproscan/main'
include { INTERPROSCAN_DATABASE } from '../../modules/local/interproscan_download'

workflow FUNCTION {
    take:
    faas // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)

    main:
    ch_versions         = Channel.empty()
    ch_interproscan_tsv = Channel.empty()
    ch_interproscan_db  = Channel.empty()

    ch_faa_for_interproscan     = faas

    if ( params.function_interproscan_db != null ) {
        ch_interproscan_db = Channel
            .fromPath( params.function_interproscan_db )
            .first() // i dont know if this is required!!!?
    } else {
        INTERPROSCAN_DATABASE ( params.function_interproscan_db_url )
        ch_versions  = ch_versions.mix( INTERPROSCAN_DATABASE.out.versions )
        ch_interproscan_db = ( INTERPROSCAN_DATABASE.out.db )
    }

    //INTERPROSCAN( ch_faa_for_interproscan, [] )
    INTERPROSCAN( ch_faa_for_interproscan, ch_interproscan_db )
    ch_interproscan_tsv = ch_interproscan_tsv.mix(INTERPROSCAN.out.tsv)

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    emit:
    versions = ch_versions
    tsv      = ch_interproscan_tsv // channel: [ val(meta), tsv ]
}
