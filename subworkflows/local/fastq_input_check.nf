/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'


workflow FASTQ_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv

    main:

    ch_versions = Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
        .tsv
        .splitCsv ( header:true, sep:'\t' )
        .map { create_fastq_channels(it) }
        .dump (tag: 'fastq_channel_before_merge_samples')
        .groupTuple(by: [0])
        .dump(tag: 'fastq_channel_after_merge_samples_grouped')
        .branch {
            meta, fastqs ->
                single: fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_reads }

    ch_versions = ch_versions.mix( SAMPLESHEET_CHECK.out.versions )

    // Merge multi-lane sample fastq for protocols except for 10x genomics, trust4 (cellranger handles multi-fastq per sample)
    if (params.library_generation_method == 'sc_10x_genomics' || params.library_generation_method == 'trust4')  {

        ch_merged_reads = ch_reads.single.mix( ch_reads.multiple )

    } else {

        CAT_FASTQ (
            ch_reads.multiple
        )
        .reads
        .mix( ch_reads.single )
        .dump (tag: 'fastq_channel_after_merge_samples')
        .set { ch_merged_reads }

        ch_versions = ch_versions.mix( CAT_FASTQ.out.versions )

    }

    emit:
    reads = ch_merged_reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
    samplesheet = SAMPLESHEET_CHECK.out.tsv // tsv metadata file
}

// Function to map
def create_fastq_channels(LinkedHashMap col) {

    def meta = [:]

    meta.id                 = col.sample_id
    meta.subject_id         = col.subject_id
    meta.species            = col.species
    meta.collapseby_group   = col."${params.collapseby}"
    meta.cloneby_group      = col."${params.cloneby}"
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
        if (!params.index_file){
            error "ERROR: --index_file was not provided but the index file path is specified in the samplesheet!"
        }
        if (!file(col.filename_I1).exists()) {
            error "ERROR: Please check input samplesheet -> Index read FastQ file does not exist!\n${col.filename_I1}"
        }
        array = [ meta, [ file(col.filename_R1), file(col.filename_R2), file(col.filename_I1) ] ]
    } else {
        array = [ meta, [ file(col.filename_R1), file(col.filename_R2) ] ]
        if (params.index_file) {
            error "ERROR: Index file path was provided but the index file path is not specified in the samplesheet!"
        }
    }
    return array
}
