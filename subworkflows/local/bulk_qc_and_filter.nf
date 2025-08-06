include { CHANGEO_CREATEGERMLINES } from '../../modules/local/changeo/changeo_creategermlines'
include { REMOVE_CHIMERIC  } from '../../modules/local/enchantr/remove_chimeric'
include { DETECT_CONTAMINATION  } from '../../modules/local/enchantr/detect_contamination'
include { COLLAPSE_DUPLICATES as COLLAPSE_DUPLICATES_WITHIN_SAMPLE } from '../../modules/local/enchantr/collapse_duplicates'
include { COLLAPSE_DUPLICATES as COLLAPSE_DUPLICATES_CROSS_SAMPLES } from '../../modules/local/enchantr/collapse_duplicates'

workflow BULK_QC_AND_FILTER {

    take:
    ch_repertoire // tuple [meta, repertoire_tab]
    ch_reference_fasta

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // Remove chimeric sequences if requested
    if (params.remove_chimeric) {

        // Create germlines (not --cloned)
        CHANGEO_CREATEGERMLINES(
            ch_repertoire,
            ch_reference_fasta.collect()
        )
        ch_logs = ch_logs.mix(CHANGEO_CREATEGERMLINES.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_CREATEGERMLINES.out.versions)

        // Remove chimera
        REMOVE_CHIMERIC(
            CHANGEO_CREATEGERMLINES.out.tab,
            ch_reference_fasta.collect()
        )
        ch_logs = ch_logs.mix(REMOVE_CHIMERIC.out.logs)
        ch_versions = ch_versions.mix(REMOVE_CHIMERIC.out.versions)
        ch_bulk_chimeric_pass = REMOVE_CHIMERIC.out.tab


    } else {
        ch_bulk_chimeric_pass = ch_repertoire
    }

    // For Bulk data, detect cross-contamination
    // This is only informative at this time
    // TODO: add a flag to specify remove suspicious sequences
    // and update file size log accordingly

    if (params.detect_contamination) {
        DETECT_CONTAMINATION(
            ch_bulk_chimeric_pass
            .map{ it -> [ it[1] ] }
            .collect()
        )
        ch_logs = ch_logs.mix(DETECT_CONTAMINATION.out.logs)
        ch_versions = ch_versions.mix(DETECT_CONTAMINATION.out.versions)
    }

    COLLAPSE_DUPLICATES_WITHIN_SAMPLE(
        ch_bulk_chimeric_pass
    )

    if (params.cross_sample_collapse == true){
            COLLAPSE_DUPLICATES_WITHIN_SAMPLE.out.tab
            .map{ it -> [   it[0].id,
                            it[0] ] }
            .set{ch_onlymeta}

        COLLAPSE_DUPLICATES_WITHIN_SAMPLE.out.tab
            .map { it -> it[1] }
            .collect()
            .map { it -> [ [id:'all_reps'], it ] }
            .set{ch_all_repertoires_within_sample}

        COLLAPSE_DUPLICATES_CROSS_SAMPLES(
            ch_all_repertoires_within_sample
        )

        COLLAPSE_DUPLICATES_CROSS_SAMPLES.out.tab
            .map { it -> it[1] }
            .flatten()
            .map { path ->
                def sample_id = path.baseName.replaceFirst("__collapse-pass", "")
                [ sample_id, path]
                }
            .set{ch_repertoire_after_2nd_collpasing_with_sampleid}

        ch_repertoire_after_2nd_collpasing_with_meta = ch_onlymeta.join(ch_repertoire_after_2nd_collpasing_with_sampleid)
                                                    .map{ it -> [ it[1], it[2] ]}

        ch_collapsed = ch_repertoire_after_2nd_collpasing_with_meta

        ch_versions = ch_versions.mix(COLLAPSE_DUPLICATES_CROSS_SAMPLES.out.versions)
        ch_logs = ch_logs.mix(COLLAPSE_DUPLICATES_CROSS_SAMPLES.out.logs)
    } else {
        ch_versions = ch_versions.mix(COLLAPSE_DUPLICATES_WITHIN_SAMPLE.out.versions)
        ch_logs = ch_logs.mix(COLLAPSE_DUPLICATES_WITHIN_SAMPLE.out.logs)
        ch_collapsed = COLLAPSE_DUPLICATES_WITHIN_SAMPLE.out.tab
    }


    emit:
    versions = ch_versions
    repertoires = ch_collapsed
    logs = ch_logs

}
