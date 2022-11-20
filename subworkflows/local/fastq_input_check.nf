/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
//include { VALIDATE_INPUT } from '../../modules/local/enchantr/validate_input'

workflow FASTQ_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .tsv
        .splitCsv ( header:true, sep:'\t' )
        .map { create_fastq_channels(it) }
        .set { ch_reads }
    // VALIDATE_INPUT(
    //     samplesheet,
    //     params.miairr,
    //     params.collapseby,
    //     params.cloneby
    // )

    // VALIDATE_INPUT.out.validated_input
    //                     .splitCsv(header: true, sep:'\t')
    //                     .map { get_meta(it) }
    //                     .set{ ch_reads }

    emit:
    reads = ch_reads // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
    samplesheet = SAMPLESHEET_CHECK.out.tsv // tsv metadata file
}

// Function to map
def create_fastq_channels(LinkedHashMap col) {

    def meta = [:]

    meta.id     = col.sample_id
    meta.filename     = col.filename
    meta.subject_id   = col.subject_id
    meta.species     = col.species
    meta.collapseby_group = col."${params.collapseby}"
    meta.cloneby_group = col."${params.cloneby}"
    meta.filetype = "fastq"
    meta.single_cell = col.single_cell
    meta.pcr_target_locus = col.pcr_target_locus
    meta.locus = col.locus

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
