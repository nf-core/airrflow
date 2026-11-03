include { CHANGEO_CREATEGERMLINES } from '../../modules/local/changeo/changeo_creategermlines'
include { REMOVE_CHIMERIC  } from '../../modules/local/enchantr/remove_chimeric'
include { DETECT_CONTAMINATION  } from '../../modules/local/enchantr/detect_contamination'
include { COLLAPSE_DUPLICATES  } from '../../modules/local/enchantr/collapse_duplicates'

workflow BULK_QC_AND_FILTER {

    take:
    ch_repertoire // tuple [meta, repertoire_tab]
    ch_reference_fasta
    remove_chimeric
    detect_contamination
    collapseby

    main:

    ch_logs = channel.empty()

    // Remove chimeric sequences if requested
    if (remove_chimeric) {

        // Create germlines (not --cloned)
        CHANGEO_CREATEGERMLINES(
            ch_repertoire,
            ch_reference_fasta.collect()
        )
        ch_logs = ch_logs.mix(CHANGEO_CREATEGERMLINES.out.logs)

        // Remove chimera
        REMOVE_CHIMERIC(
            CHANGEO_CREATEGERMLINES.out.tab,
            ch_reference_fasta.collect()
        )
        ch_logs = ch_logs.mix(REMOVE_CHIMERIC.out.logs)
        ch_bulk_chimeric_pass = REMOVE_CHIMERIC.out.tab


    } else {
        ch_bulk_chimeric_pass = ch_repertoire
    }

    // For Bulk data, detect cross-contamination
    // This is only informative at this time
    // TODO: add a flag to specify remove suspicious sequences
    // and update file size log accordingly

    if (detect_contamination) {
        DETECT_CONTAMINATION(
            ch_bulk_chimeric_pass
            .map{ it -> [ it[1] ] }
            .collect()
        )
        ch_logs = ch_logs.mix(DETECT_CONTAMINATION.out.logs)
    }

    COLLAPSE_DUPLICATES(
        ch_bulk_chimeric_pass,
        collapseby
    )

    ch_logs = ch_logs.mix(COLLAPSE_DUPLICATES.out.logs)

    emit:
    repertoires = COLLAPSE_DUPLICATES.out.tab
    logs = ch_logs

}
