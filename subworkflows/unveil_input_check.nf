/*
 * Check input samplesheet and get channels
 */

params.options = [:]

include {
    VALIDATE_INPUT
    } from '../modules/local/unveil/validate_input' addParams( options: params.options )

workflow UNVEIL_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    miairr      
    collapseby
    cloneby

    main:
    VALIDATE_INPUT ( samplesheet, miairr, collapseby, cloneby).validated_input
            .splitCsv(header: true, sep:'\t')
            .map { get_meta(it) }
            .branch { it ->
                fasta: it.filename =~ /[fasta|fa]$/
                tsv:   it.filename =~ /tsv$/
            }
            .set{ch_metadata}

    emit:
    ch_fasta = ch_metadata.fasta
    ch_tsv = ch_metadata.tsv
    //ch_metadata
}


// Function to map 
def get_meta (LinkedHashMap col) {

    def meta = [:]

    meta.input_id     = col.input_id
    meta.filename     = file(col.filename, checkifExists: true)
    meta.subject_id   = col.subject_id
    meta.organism     = col.organism
    meta.collapseby_group = col.collapseby_group
    meta.collapseby_size  = col.collapseby_size
    meta.cloneby_group = col.cloneby_group
    meta.cloneby_size = col.cloneby_size
    meta.filetype = col.filetype

    return meta    
}