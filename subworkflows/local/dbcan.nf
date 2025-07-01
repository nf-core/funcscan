/*
    Run rundbcan screening tools
*/

include { RUNDBCAN_DATABASE         } from '../../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_CAZYMEANNOTATION } from '../../modules/nf-core/rundbcan/cazymeannotation/main'
include { RUNDBCAN_EASYCGC          } from '../../modules/nf-core/rundbcan/easycgc/main'
include { RUNDBCAN_EASYSUBSTRATE    } from '../../modules/nf-core/rundbcan/easysubstrate/main'


workflow DBCAN {

    take:
    faas    // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)
    gffs    // tuple val(meta), path(ANNOTATION_ANNOTATION_TOOL.out.gff)

    main:

    ch_versions = Channel.empty()

    // When adding new tool that requires FAA, make sure to update conditions
    // in funcscan.nf around annotation and AMP subworkflow execution
    // to ensure annotation is executed!
    ch_faas_for_rundbcan = faas
    ch_gffs_for_rundbcan = gffs

    // RUN DBCAN
    RUNDBCAN_DATABASE ()
    ch_versions = ch_versions.mix(RUNDBCAN_DATABASE.out.versions)

    // RUN CAZyme Annotation
    RUNDBCAN_CAZYMEANNOTATION (
        ch_faas_for_rundbcan,
        RUNDBCAN_DATABASE.out.dbcan_db
    )
    ch_versions = ch_versions.mix(RUNDBCAN_CAZYMEANNOTATION.out.versions)

    // Prepare input for DBCAN CGC and SUBSTRATE
    ch_input_for_dbcan = ch_faas_for_rundbcan
        .join(ch_gffs_for_rundbcan)
        .multiMap { meta, faa, gff ->
            faa: [meta, faa]
            gff: [meta, gff, 'prodigal']
        }

    // Run DBCAN CGC Annotation
    if ( !params.dbcan_skip_cgc ) {
        RUNDBCAN_EASYCGC (
            ch_input_for_dbcan.faa,
            ch_input_for_dbcan.gff,
            RUNDBCAN_DATABASE.out.dbcan_db
        )
        ch_versions = ch_versions.mix(RUNDBCAN_EASYCGC.out.versions)
    }

    // Run DBCAN Substrate Annotation
    if ( !params.dbcan_skip_substrate ) {
        RUNDBCAN_EASYSUBSTRATE (
            ch_input_for_dbcan.faa,
            ch_input_for_dbcan.gff,
            RUNDBCAN_DATABASE.out.dbcan_db
        )
        ch_versions = ch_versions.mix(RUNDBCAN_EASYSUBSTRATE.out.versions)
    }

    emit:
    versions = ch_versions
}
