/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { SAMPLESHEET_CHECK } from '../process/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv
    
    main:
    // TODO: avoiding checking samplesheet for now, add samplesheet check later.
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
    meta.population   = col.Population

    def array = []
    if (!file(col.R1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${col.fastq_1}"
    }
    if (col.I1) {
        if (!file(col.I1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Index read FastQ file does not exist!\n${col.fastq_2}"
        }
        array = [ meta, [ file(col.R1), file(col.R2), file(col.I1) ] ]
    } else {
        if (!file(col.R2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${col.fastq_2}"
        }
        array = [ meta, [ file(col.R1), file(col.R2) ] ]
    }
    return array    
}