/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow FASTQ_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .tsv
        .splitCsv ( header:true, sep:'\t' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    samplesheet = SAMPLESHEET_CHECK.out.tsv // tsv metadata file
}

// Function to map
def create_fastq_channels(LinkedHashMap col) {
    def meta = [:]
    meta.id           = col.sample_id
    meta.subject      = col.subject_id
    meta.locus        = col.pcr_target_locus
    meta.species      = col.species
    meta.single_cell  = 'false'

    def array = []
    if (!file(col.filename_R1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${col.filename_R1}"
    }
    if (!file(col.filename_R2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${col.filename_R2}"
    }
    if (col.filename_I1) {
        if (!file(col.filename_I1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Index read FastQ file does not exist!\n${col.filename_I1}"
        }
        array = [ meta, [ file(col.filename_R1), file(col.filename_R2), file(col.filename_I1) ] ]
    } else {

        array = [ meta, [ file(col.filename_R1), file(col.filename_R2) ] ]
    }
    return array
}
