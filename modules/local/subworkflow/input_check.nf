/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include {
    SAMPLESHEET_CHECK;
    get_samplesheet_paths } from '../process/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv
    
    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:'\t' )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to map 
def get_samplesheet_paths(LinkedHashMap col) {
    def meta = [:]
    meta.id           = col.ID
    meta.source       = col.Source
    meta.treatment    = col.Treatment
    meta.time         = col.Extraction_time
    meta.single_end   = collect.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array    
}