include { SINGLE_CELL_QC  } from '../../modules/local/enchantr/single_cell_qc'

workflow SINGLE_CELL_QC_AND_FILTERING {
    take:
    repertoires // tuple [meta, repertoire_tab]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    repertoires
            .map{ it -> [   it[0].id,
                            it[0] ] }
            .set{ch_onlymeta}

    repertoires
            .map { it -> it[1]}
            .collect()
            .set{ch_repertoire_allsc}

    SINGLE_CELL_QC(
        ch_repertoire_allsc
    )

    SINGLE_CELL_QC.out.tab
                .flatten()
                .map { it -> [ "${it.baseName}".replaceFirst("__scqc-pass", ""), it ] }
                .set{ch_repertoire_after_scqc_with_sampleid}

    ch_logs = ch_logs.mix(SINGLE_CELL_QC.out.logs)
    ch_versions = ch_versions.mix(SINGLE_CELL_QC.out.versions.ifEmpty(null))

    ch_repertoire_after_scqc_withmeta = ch_onlymeta.join(ch_repertoire_after_scqc_with_sampleid)
                                                    .map{ it -> [ it[1], it[2] ]}

    emit:
    versions = ch_versions
    repertoires = ch_repertoire_after_scqc_withmeta
    logs = ch_logs
}
