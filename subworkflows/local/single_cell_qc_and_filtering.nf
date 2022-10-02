include { SINGLE_CELL_QC  } from '../../modules/local/enchantr/single_cell_qc'

workflow SINGLE_CELL_QC_AND_FILTERING {
    take:
    repertoires // tuple [meta, repertoire_tab]

    main:
    ch_versions = Channel.empty()

    SINGLE_CELL_QC(
        repertoires
        .map{ it -> [ it[1] ] }
        .collect()
    )
    // ch_file_sizes = ch_file_sizes.mix(SINGLE_CELL_QC.out.logs)
    ch_versions = ch_versions.mix(SINGLE_CELL_QC.out.versions.ifEmpty(null))

    emit:
    versions = ch_versions
    repertoires = SINGLE_CELL_QC.out.tab
}
