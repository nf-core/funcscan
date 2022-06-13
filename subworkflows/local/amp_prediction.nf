//
// Check input samplesheet and get read channels
//

include { MACREL_CONTIGS } from '../../modules/nf-core/macrel/contigs/main'
include { AMPIR_PEPTIDES } from '../../modules/foo'
include { MERGE_PREDICTIONS } from '../../modules/local/merge_predictions'

workflow AMP_PREDICTION {
    take:
    proteins  // tuple val(meta), path('*.faa')

    main:
    MACREL_PEPTIDES(proteins)
    AMPIR_PEPTIDES(proteins)
    // AMPLIFY

    ch_predictions = MACREL_CONTIGS.out.amp_prediction.join(
        AMPIR_PEPTIDES.out.amp_prediction
    )
    // .join(...)
    MERGE_PREDICTIONS(ch_predictions)

    emit:
    predictions = MERGE_PREDICTIONS.out.                                   // channel: [ val(meta), [ fasta ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
