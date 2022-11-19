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

    ch_repertoire.map{ it -> [ it[0]."${params.cloneby}", it[0], it[1] ] }
                .dump(tag:'cloneby')

    DEFINE_CLONES(
        ch_repertoire,
        clone_threshold,
        ch_imgt
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
