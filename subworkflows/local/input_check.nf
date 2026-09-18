/*
 * Check the input samplesheet, validate it and stage the input files.
 *
 * This subworkflow is run at the very beginning of the pipeline, before the
 * germline references are fetched, so that the sample metadata (species,
 * locus, subject_id) is available to every downstream subworkflow.
 */

include { SAMPLESHEET_CHECK                                } from '../../modules/local/samplesheet_check'
include { SAMPLESHEET_CHECK as SAMPLESHEET_CHECK_ASSEMBLED } from '../../modules/local/samplesheet_check'
include { VALIDATE_INPUT                                   } from '../../modules/local/enchantr/validate_input'
include { CAT_FASTQ                                        } from '../../modules/nf-core/cat/fastq/main'
include { RENAME_FILE as RENAME_FILE_FASTA                 } from '../../modules/local/rename_file'
include { RENAME_FILE as RENAME_FILE_TSV                   } from '../../modules/local/rename_file'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv
    mode
    library_generation_method
    miairr
    collapseby
    cloneby
    reassign
    index_file

    main:

    ch_versions = channel.empty()
    ch_logs     = channel.empty()

    ch_reads = channel.empty()
    ch_fasta = channel.empty()
    ch_tsv   = channel.empty()

    if ( mode == "fastq" ) {

        SAMPLESHEET_CHECK ( samplesheet )
            .tsv
            .splitCsv ( header:true, sep:'\t' )
            .map { create_fastq_channels(it, collapseby, cloneby, index_file) }
            .groupTuple(by: [0])
            .branch {
                meta, fastqs ->
                    single: fastqs.size() == 1
                        return [ meta, fastqs.flatten() ]
                    multiple: fastqs.size() > 1
                        return [ meta, fastqs.flatten() ]
            }
            .set { ch_split_reads }

        // Merge multi-lane sample fastq for protocols except for 10x genomics, trust4 (cellranger handles multi-fastq per sample)
        if (library_generation_method == 'sc_10x_genomics' || library_generation_method == 'trust4')  {

            ch_reads = ch_split_reads.single.mix( ch_split_reads.multiple )

        } else {

            CAT_FASTQ (
                ch_split_reads.multiple
            )
            .reads
            .mix( ch_split_reads.single )
            .dump (tag: 'fastq_channel_after_merge_samples')
            .set { ch_reads }

            ch_versions = ch_versions.mix( CAT_FASTQ.out.versions )
        }

        ch_validated_samplesheet = SAMPLESHEET_CHECK.out.tsv

    } else if ( mode == "assembled" ) {

        SAMPLESHEET_CHECK_ASSEMBLED ( samplesheet )
        VALIDATE_INPUT (
            samplesheet,
            miairr,
            collapseby,
            cloneby,
            reassign
        )
        ch_validated_samplesheet = VALIDATE_INPUT.out.validated_input

        ch_validated_samplesheet
            .splitCsv(header: true, sep:'\t')
            .map { get_meta(it) }
                .branch { it ->
                    fasta: it[0].filename =~ /[fasta|fa]$/
                    tsv:   it[0].filename =~ /tsv$/
                }
                .set{ ch_metadata }

        RENAME_FILE_FASTA( ch_metadata.fasta )
        ch_fasta = RENAME_FILE_FASTA.out.file
        ch_logs = ch_logs.mix(RENAME_FILE_FASTA.out.logs)

        RENAME_FILE_TSV( ch_metadata.tsv )
        ch_tsv = RENAME_FILE_TSV.out.file
        ch_logs = ch_logs.mix(RENAME_FILE_TSV.out.logs)

    } else {
        error "Mode parameter value not valid."
    }

    emit:
    reads = ch_reads // channel: [ val(meta), [ reads ] ], fastq mode only
    fasta = ch_fasta // channel: [ val(meta), fasta ], assembled mode only
    tsv = ch_tsv // channel: [ val(meta), tsv ], assembled mode only
    validated_samplesheet = ch_validated_samplesheet // tsv metadata file
    logs = ch_logs // channel: [ logs ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to map the raw (fastq) samplesheet
def create_fastq_channels(LinkedHashMap col, collapseby, cloneby, index_file) {

    def meta = [:]

    meta.id                 = col.sample_id
    meta.sample_id          = col.sample_id
    meta.subject_id         = col.subject_id
    meta.species            = col.species
    meta.collapseby_group   = col[collapseby]
    meta.cloneby_group      = col[cloneby]
    meta.filetype           = "fastq"
    meta.single_cell        = col.single_cell.toLowerCase()
    meta.locus              = col.pcr_target_locus
    meta.single_end         = false

    def array = []
    if (!file(col.filename_R1).exists()) {
        error "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${col.filename_R1}"
    }
    if (!file(col.filename_R2).exists()) {
        error "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${col.filename_R2}"
    }
    if (col.filename_I1) {
        if (!index_file){
            error "ERROR: --index_file was not provided but the index file path is specified in the samplesheet!"
        }
        if (!file(col.filename_I1).exists()) {
            error "ERROR: Please check input samplesheet -> Index read FastQ file does not exist!\n${col.filename_I1}"
        }
        array = [ meta, [ file(col.filename_R1), file(col.filename_R2), file(col.filename_I1) ] ]
    } else {
        array = [ meta, [ file(col.filename_R1), file(col.filename_R2) ] ]
        if (index_file) {
            error "ERROR: Index file path was provided but the index file path is not specified in the samplesheet!"
        }
    }
    return array
}

// Function to map the validated (assembled) samplesheet
def get_meta (LinkedHashMap col) {

    def meta = [:]

    meta.id     = col.sample_id
    meta.filename     = col.filename
    meta.sample_id = col.sample_id
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
