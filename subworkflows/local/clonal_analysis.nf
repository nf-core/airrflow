include { FIND_THRESHOLD  } from '../../modules/local/enchantr/find_threshold'
include { DEFINE_CLONES  } from '../../modules/local/enchantr/define_clones'
include { DOWSER_LINEAGES } from '../../modules/local/enchantr/dowser_lineages'

workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire
    ch_imgt
    ch_logo

    main:
    ch_versions = Channel.empty()


    if (params.clonal_threshold == "auto") {

        ch_find_threshold = ch_repertoire.map{ it -> it[1] }
                                        .collect()
                                        .dump(tag:'find_threshold')

        FIND_THRESHOLD (
            ch_find_threshold,
            ch_logo
        )
        ch_threshold = FIND_THRESHOLD.out.mean_threshold
        ch_versions = ch_versions.mix(FIND_THRESHOLD.out.versions)

        clone_threshold = ch_threshold
            .splitText( limit:1 ) { it.trim().toString() }
            .dump(tag: 'clone_threshold')
            .filter { it != 'NA'}
            .filter { it != 'NaN' }
            .dump(tag: "threshold")
            .ifEmpty { exit 1, "Automatic clone_threshold is 'NA'. Consider setting params.threshold manually."}

    } else {
        clone_threshold = params.clonal_threshold
    }

    ch_repertoire.map{ it -> [ it[0]."${params.cloneby}",
                                it[0].id,
                                it[0].filename,
                                it[0].subject_id,
                                it[0].species,
                                it[0].filetype,
                                it[0].single_cell,
                                it[0].pcr_target_locus,
                                it[0].locus,
                                it[1] ] }
                .groupTuple()
                .dump(tag:'cloneby')
                .map{ get_meta_tabs(it) }
                .dump(tag:'cloneby_after_map')
                .set{ ch_define_clones }

    DEFINE_CLONES(
        ch_define_clones,
        clone_threshold.collect(),
        ch_imgt.collect()
    )
    ch_versions = ch_versions.mix(DEFINE_CLONES.out.versions)

    if (!params.skip_lineage){
        DOWSER_LINEAGES(
            DEFINE_CLONES.out.tab
                .flatten()
        )
        ch_versions = ch_versions.mix(DOWSER_LINEAGES.out.versions)
    }

    emit:
    repertoire = DEFINE_CLONES.out.tab
    versions = ch_versions
}

// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.cloneby            = [arr[0]].unique().join("")
    meta.sample_ids         = arr[1]
    meta.filename           = arr[2]
    meta.subject_id         = arr[3]
    meta.species            = arr[4]
    meta.filetype           = arr[5]
    meta.single_cell        = arr[6].unique().join("")
    meta.pcr_target_locus   = arr[7].unique().join("")
    meta.locus              = arr[8]

    def array = []

        array = [ meta, arr[9].flatten() ]

    return array
}
