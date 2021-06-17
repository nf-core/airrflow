/*
 * Check input samplesheet and get channels
 */

params.options = [:]

include {
    VALIDATE_INPUT
    } from '../modules/local/airrflow/validate_input' addParams( options: params.options )

workflow AIRRFLOW_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    miairr      
    collapseby
    cloneby

    main:
    VALIDATE_INPUT ( samplesheet, miairr, collapseby, cloneby).validated_input
              .splitCsv(header: true, sep:'\t')
              .map { row -> tuple(
                                 file("${row.filename}", checkifExists: true), 
                                 "${row.subject_id}",
                                 "${row.organism}",
                                 "${row.collapseby_group}",
                                 "${row.collapseby_size}",
                                 "${row.cloneby_group}",
                                 "${row.cloneby_size}",
                                 "${row.filetype}",
                                 "${row.input_id}"
                                 )
              }
             .branch { filename, subject_id, organism, collapseby_group, collapseby_size, cloneby_group, cloneby_size, filetype, input_id ->
                fasta: "${filename}" =~ /[fasta|fa]$/
                tsv:   "${filename}" =~ /tsv$/
              }
              .set{ch_metadata}

    emit:
    ch_fasta = ch_metadata.fasta
    ch_tsv = ch_metadata.tsv
}
