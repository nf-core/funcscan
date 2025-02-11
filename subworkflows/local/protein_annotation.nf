/*
    RUN FUNCTIONAL CLASSIFICATION
*/

include { INTERPROSCAN_DATABASE } from '../../modules/local/interproscan_download'
include { INTERPROSCAN          } from '../../modules/nf-core/interproscan/main'

workflow PROTEIN_ANNOTATION {
    take:
    faas // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)

    main:
    ch_versions                  = Channel.empty()
    ch_interproscan_tsv          = Channel.empty()
    ch_interproscan_db           = Channel.empty()
    ch_interproscan_tsv_modified = Channel.empty()

    ch_faa_for_interproscan = faas

    if ( params.protein_annotation_tool == 'InterProScan') {

        if ( params.protein_annotation_interproscan_db != null ) {
            ch_interproscan_db = Channel
                .fromPath( params.protein_annotation_interproscan_db )
                .first()
        } else {
            INTERPROSCAN_DATABASE ( params.protein_annotation_interproscan_db_url )
            ch_versions = ch_versions.mix( INTERPROSCAN_DATABASE.out.versions )
            ch_interproscan_db = ( INTERPROSCAN_DATABASE.out.db )
        }

        INTERPROSCAN( ch_faa_for_interproscan, ch_interproscan_db )
        ch_versions = ch_versions.mix( INTERPROSCAN.out.versions )
        ch_interproscan_tsv = ch_interproscan_tsv.mix( INTERPROSCAN.out.tsv )

        // Current INTERPROSCAN version 5.59_91.0 only includes 13 columns and not 15 which ampcombi expects, so we added them here
        ch_interproscan_tsv_modified = INTERPROSCAN.out.tsv
            .map { meta, tsv_path ->
                def modified_tsv_path = "${workflow.workDir}/tmp/${meta.id}_interproscan.faa.tsv"

                def modified_tsv_content = new File(tsv_path.toString())
                    .readLines()
                    .collect { line -> (line.split('\t') + ['NA', 'NA']).join('\t') }

                new File(modified_tsv_path).text = modified_tsv_content.join('\n')
                [meta, file(modified_tsv_path)]
            }

        ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)
    }

    emit:
    versions = ch_versions
    tsv      = ch_interproscan_tsv_modified // channel: [ val(meta), tsv ]
}
