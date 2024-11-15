include { CHANGEO_ASSIGNGENES } from '../../modules/local/changeo/changeo_assigngenes'
include { CHANGEO_MAKEDB } from '../../modules/local/changeo/changeo_makedb'
include { CHANGEO_PARSEDB_SPLIT } from '../../modules/local/changeo/changeo_parsedb_split'
// reveal
include { FILTER_QUALITY  } from '../../modules/local/reveal/filter_quality'
include { FILTER_JUNCTION_MOD3  } from '../../modules/local/reveal/filter_junction_mod3'
include { ADD_META_TO_TAB  } from '../../modules/local/reveal/add_meta_to_tab'


workflow VDJ_ANNOTATION {

    take:
    ch_fasta // [meta, fasta]
    ch_validated_samplesheet
    ch_igblast
    ch_reference_fasta

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    CHANGEO_ASSIGNGENES (
        ch_fasta,
        ch_igblast.collect()
    )

    ch_logs = ch_logs.mix(CHANGEO_ASSIGNGENES.out.logs)
    ch_versions = ch_versions.mix(CHANGEO_ASSIGNGENES.out.versions)

    CHANGEO_MAKEDB (
        CHANGEO_ASSIGNGENES.out.fasta,
        CHANGEO_ASSIGNGENES.out.blast,
        ch_reference_fasta.collect()
    )
    ch_logs = ch_logs.mix(CHANGEO_MAKEDB.out.logs)
    ch_versions = ch_versions.mix(CHANGEO_MAKEDB.out.versions)

    ch_assigned_tab = CHANGEO_MAKEDB.out.tab
    ch_assignment_logs = CHANGEO_MAKEDB.out.logs

    if (!params.skip_alignment_filter){
        // Apply quality filters:
        // - locus should match v_call chain
        // - seq alignment min length informative positions 200
        // - max 10% N nucleotides
        FILTER_QUALITY(
            ch_assigned_tab
        )
        ch_for_parsedb_split = FILTER_QUALITY.out.tab
        ch_logs = ch_logs.mix(FILTER_QUALITY.out.logs)
        ch_versions = ch_versions.mix(FILTER_QUALITY.out.versions)
    } else {
        ch_for_parsedb_split = ch_assigned_tab
    }

    if (params.productive_only) {
        CHANGEO_PARSEDB_SPLIT (
            ch_for_parsedb_split
        )
        ch_logs = ch_logs.mix(CHANGEO_PARSEDB_SPLIT.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_PARSEDB_SPLIT.out.versions)

        // Apply filter: junction length multiple of 3
        FILTER_JUNCTION_MOD3(
            CHANGEO_PARSEDB_SPLIT.out.tab
        )
        ch_logs = ch_logs.mix(FILTER_JUNCTION_MOD3.out.logs)
        ch_versions = ch_versions.mix(FILTER_JUNCTION_MOD3.out.versions)
        ch_repertoire = FILTER_JUNCTION_MOD3.out.tab

    } else {
        ch_repertoire = FILTER_QUALITY.out.tab
    }

    ADD_META_TO_TAB(
        ch_repertoire,
        ch_validated_samplesheet
    )
    ch_logs = ch_logs.mix(ADD_META_TO_TAB.out.logs)
    ch_versions = ch_versions.mix(ADD_META_TO_TAB.out.versions)


    emit:
    versions = ch_versions
    repertoire = ADD_META_TO_TAB.out.tab
    reference_fasta = ch_reference_fasta
    reference_igblast = ch_igblast
    changeo_makedb_logs = ch_assignment_logs
    logs = ch_logs

}
