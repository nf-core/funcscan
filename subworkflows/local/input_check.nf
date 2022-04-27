//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    parsed_samplesheet = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channels(it) }

    emit:
    contigs  = parsed_samplesheet             // channel: [ val(meta), [ fasta ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fasta ] ]
def create_input_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    if (!file(row.fasta).exists()) {
        exit 1, "[funcscan] error: please check input samplesheet. FASTA file does not exist for: \n${row.fasta}"
    } else if ( row.faa != '' && !file(row.faa).exists() ) {
        exit 1, "[funcscan] error: please check input samplesheet. FAA file does not exist for: \n${row.faa}"
    } else {
        def faafile = row.faa != '' ? file(row.faa) : ''
        array = [ meta, file(row.fasta), faafile ]
    }

    return array
}
