/*
 * Check input samplesheet and get channels
 */

include { VALIDATE_INPUT } from '../../modules/local/enchantr/validate_input'

workflow ASSEMBLED_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    miairr
    collapseby
    cloneby

    main:
    // TODO: validate input should check that sample_ids are unique

    VALIDATE_INPUT ( samplesheet, miairr, collapseby, cloneby ) //removed reassign
    ch_validated_input = VALIDATE_INPUT.out.validated_input
    ch_validated_input
        .splitCsv(header: true, sep:'\t')
        .map { get_meta(it) }
            .branch { it ->
                fasta: it[0].filename =~ /[fasta|fa]$/
                tsv:   it[0].filename =~ /tsv$/
            }
            .set{ ch_metadata }

    emit:
    ch_fasta = ch_metadata.fasta
    ch_tsv = ch_metadata.tsv
    validated_input = ch_validated_input
    versions = VALIDATE_INPUT.out.versions
}

// Function to map
def get_meta (LinkedHashMap col) {

    def meta = [:]

    meta.id     = col.sample_id
    meta.filename     = col.filename
    meta.subject_id   = col.subject_id
    meta.species     = col.species
    meta.collapseby_group = col.collapseby_group
    meta.collapseby_size  = col.collapseby_size
    meta.cloneby_group = col.cloneby_group
    meta.cloneby_size = col.cloneby_size
    meta.filetype = col.filetype
    meta.single_cell = col.single_cell
    meta.pcr_target_locus = col.pcr_target_locus
    meta.locus = col.locus

    if (!file(col.filename).exists()) {
        error "ERROR: Please check input samplesheet: filename does not exist!\n${col.filename}"
    }

    return  [ meta, file(col.filename) ]
}
