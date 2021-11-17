//
// Check input samplesheet and get read channels
//

params.options = [:]

include { GUNZIP } from '../../modules/nf-core/modules/gunzip/main' addParams( options: params.options )

workflow GUNZIP_WITH_META {
    take:
    input // input tuple [ val(meta), path(fasta.gz) ]

    main:
    input.dump(tag: "hello")

    GUNZIP ( input )
        .set { output }

    emit:
    output = [ input[1], GUNZIP.out.gunzip ]     // channel: [ val(meta), [ fasta ] ]
    versions = GUNZIP.out.versions // channel: [ versions.yml ]
}
