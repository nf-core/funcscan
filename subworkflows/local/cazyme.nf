/*
    Run rundbcan screening tools
*/

include { RUNDBCAN_DATABASE         } from '../../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_CAZYMEANNOTATION } from '../../modules/nf-core/rundbcan/cazymeannotation/main'
include { RUNDBCAN_EASYCGC          } from '../../modules/nf-core/rundbcan/easycgc/main'
include { RUNDBCAN_EASYSUBSTRATE    } from '../../modules/nf-core/rundbcan/easysubstrate/main'


workflow CAZYME {
    take:
    faas // tuple val(meta), path(PROKKA/PRODIGAL.out.faa)
    gffs // tuple val(meta), path(ANNOTATION_ANNOTATION_TOOL.out.gff)

    main:
    // When adding new tool that requires FAA, make sure to update conditions
    // in funcscan.nf around annotation and dbCAN subworkflow execution
    // to ensure annotation is executed!
    ch_faas_for_rundbcan = faas
    ch_gffs_for_rundbcan = gffs

    // Prepare channel for database
    if (!params.cazyme_skip_dbcan && params.cazyme_dbcan_db) {
        ch_dbcan_db = Channel
            .fromPath(params.cazyme_dbcan_db, checkIfExists: true)
            .first()
    }
    else if (!params.cazyme_skip_dbcan && !params.cazyme_dbcan_db) {
        // Download dbCAN database
        RUNDBCAN_DATABASE()
        ch_dbcan_db = RUNDBCAN_DATABASE.out.dbcan_db
    }

    if (!params.cazyme_skip_dbcan) {
        // CAZyme annotation
        RUNDBCAN_CAZYMEANNOTATION(ch_faas_for_rundbcan, ch_dbcan_db)

        // Prepare input for dbCAN CGC and substrate annotation
        if (!params.dbcan_skip_cgc || !params.dbcan_skip_substrate) {
            ch_input_for_dbcan = ch_faas_for_rundbcan
                .join(ch_gffs_for_rundbcan)
                .filter { meta, faa, gff ->
                    if (!gff || !meta.gff_type) {
                        log.warn("Skipping sample: ${meta.id ?: 'unknown'} for dbcan CGC and substrate annotation due to empty gff or gff_type")
                        return false
                    }
                    return true
                }
                .multiMap { meta, faa, gff ->
                    faa: [meta, faa]
                    gff: [meta, gff, meta.gff_type]
                }

            // CGC annotation
            if (!params.dbcan_skip_cgc) {
                RUNDBCAN_EASYCGC(
                    ch_input_for_dbcan.faa,
                    ch_input_for_dbcan.gff,
                    ch_dbcan_db,
                )
            }


            // substrate annotation
            if (!params.dbcan_skip_substrate) {
                RUNDBCAN_EASYSUBSTRATE(
                    ch_input_for_dbcan.faa,
                    ch_input_for_dbcan.gff,
                    ch_dbcan_db,
                )
            }
        }
    }
}
