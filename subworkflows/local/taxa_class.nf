/*
    TAXONOMIC CLASSIFICATION
*/

include { MMSEQS_CREATEDB  } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_DATABASES } from '../../modules/nf-core/mmseqs/databases/main'
include { MMSEQS_TAXONOMY  } from '../../modules/nf-core/mmseqs/taxonomy/main'
include { MMSEQS_CREATETSV } from '../../modules/nf-core/mmseqs/createtsv/main'

workflow TAXA_CLASS {
    take:
    contigs // tuple val(meta), path(contigs)

    main:
    ch_versions = Channel.empty()
    ch_mmseqs_db = Channel.empty()
    ch_taxonomy_querydb = Channel.empty()
    ch_taxonomy_querydb_taxdb = Channel.empty()
    ch_taxonomy_tsv = Channel.empty()

    if (params.taxa_classification_tool == 'mmseqs2') {

        // Download the ref db if not supplied by user
        // MMSEQS_DATABASE
        if (params.taxa_classification_mmseqs_db != null) {
            ch_mmseqs_db = Channel
                .fromPath(params.taxa_classification_mmseqs_db, checkIfExists: true)
                .first()
        }
        else {
            MMSEQS_DATABASES(params.taxa_classification_mmseqs_db_id)
            ch_versions = ch_versions.mix(MMSEQS_DATABASES.out.versions)
            ch_mmseqs_db = MMSEQS_DATABASES.out.database
        }

        // Create db for query contigs, assign taxonomy and convert to table format
        // MMSEQS_CREATEDB
        MMSEQS_CREATEDB(contigs)
        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)

        // MMSEQS_TAXONOMY
        MMSEQS_TAXONOMY(MMSEQS_CREATEDB.out.db, ch_mmseqs_db)
        ch_versions = ch_versions.mix(MMSEQS_TAXONOMY.out.versions)
        ch_taxonomy_querydb_taxdb = MMSEQS_TAXONOMY.out.db_taxonomy

        // Join together to ensure in sync
        ch_taxonomy_input_for_createtsv = MMSEQS_CREATEDB.out.db
            .join(MMSEQS_TAXONOMY.out.db_taxonomy)
            .multiMap { meta, db, db_taxonomy ->
                db: [meta, db]
                taxdb: [meta, db_taxonomy]
            }

        // MMSEQS_CREATETSV
        MMSEQS_CREATETSV(ch_taxonomy_input_for_createtsv.taxdb, [[:], []], ch_taxonomy_input_for_createtsv.db)
        ch_versions = ch_versions.mix(MMSEQS_CREATETSV.out.versions)
        ch_taxonomy_tsv = MMSEQS_CREATETSV.out.tsv
    }

    emit:
    versions        = ch_versions
    sample_taxonomy = ch_taxonomy_tsv // channel: [ val(meta), tsv ]
}
