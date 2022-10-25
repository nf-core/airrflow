include { CHANGEO_CREATEGERMLINES } from '../../modules/local/changeo/changeo_creategermlines'
include { REMOVE_CHIMERIC  } from '../../modules/local/enchantr/remove_chimeric'
include { DETECT_CONTAMINATION  } from '../../modules/local/enchantr/detect_contamination'
include { COLLAPSE_DUPLICATES  } from '../../modules/local/enchantr/collapse_duplicates'

workflow BULK_QC_AND_FILTER {

    take:
    ch_repertoire // tuple [meta, repertoire_tab]
    ch_imgt

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // Remove chimeric sequences if requested
    if (params.remove_chimeric) {

        // Create germlines (not --cloned)
        CHANGEO_CREATEGERMLINES(
            ch_repertoire,
            ch_imgt.collect()
        )
        ch_logs = ch_logs.mix(CHANGEO_CREATEGERMLINES.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_CREATEGERMLINES.out.versions.ifEmpty(null))

        // Remove chimera
        REMOVE_CHIMERIC(
            CHANGEO_CREATEGERMLINES.out.tab,
            ch_imgt.collect()
        )
        ch_logs = ch_logs.mix(REMOVE_CHIMERIC.out.logs)
        ch_versions = ch_versions.mix(REMOVE_CHIMERIC.out.versions.ifEmpty(null))
        ch_bulk_chimeric_pass = REMOVE_CHIMERIC.out.tab


    } else {
        ch_bulk_chimeric_pass = ch_repertoire
    }

    // For Bulk data, detect cross-contamination
    // This is only informative at this time
    // TODO: add a flag to specify remove suspicious sequences
    // and update file size log accordingly

    DETECT_CONTAMINATION(
        ch_bulk_chimeric_pass
            .map{ it -> [ it[1] ] }
            .collect()
    )
    // TODO file size
    ch_versions = ch_versions.mix(DETECT_CONTAMINATION.out.versions.ifEmpty(null))

    ch_for_collapse = ch_bulk_chimeric_pass
            .dump()

    COLLAPSE_DUPLICATES(
        //ch_bulk_chimeric_pass
        //    .map{ it -> [ it[1] ] }
        //    .collect()
        ch_for_collapse
    )

    ch_versions = ch_versions.mix(COLLAPSE_DUPLICATES.out.versions.ifEmpty(null))
    // TODO file size

    emit:
    versions = ch_versions
    repertoires = COLLAPSE_DUPLICATES.out.tab
    logs = ch_logs

}
