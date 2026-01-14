/*
 * Check input samplesheet and get channels
 */

include { VALIDATE_INPUT } from '../../modules/local/enchantr/validate_input'
include { SAMPLESHEET_CHECK as SAMPLESHEET_CHECK_ASSEMBLED } from '../../modules/local/samplesheet_check'
include { RENAME_FILE as RENAME_FILE_FASTA } from '../../modules/local/rename_file'
include { RENAME_FILE as RENAME_FILE_TSV } from '../../modules/local/rename_file'

workflow ASSEMBLED_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    miairr
    collapseby
    cloneby

    main:
    ch_logs = Channel.empty()

    SAMPLESHEET_CHECK_ASSEMBLED ( samplesheet )
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

    RENAME_FILE_FASTA( ch_metadata.fasta )
    ch_unique_fasta = RENAME_FILE_FASTA.out.file
    ch_logs = ch_logs.mix(RENAME_FILE_FASTA.out.logs)

    RENAME_FILE_TSV( ch_metadata.tsv )
    ch_unique_tsv = RENAME_FILE_TSV.out.file
    ch_logs = ch_logs.mix(RENAME_FILE_TSV.out.logs)

    emit:
    ch_fasta = ch_unique_fasta
    ch_tsv = ch_unique_tsv
    validated_input = ch_validated_input
    versions = VALIDATE_INPUT.out.versions
    logs = ch_logs
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
