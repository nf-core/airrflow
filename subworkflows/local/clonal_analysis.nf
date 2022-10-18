include { FIND_THRESHOLD  } from '../../modules/local/enchantr/find_threshold'
include { DEFINE_CLONES  } from '../../modules/local/enchantr/define_clones'

workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire
    ch_imgt

    main:
    ch_versions = Channel.empty()

    if (params.clonal_threshold == "auto") {
        FIND_THRESHOLD (
            ch_repertoire
        )
        ch_threshold = FIND_THRESHOLD.out.mean_threshold

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
    //TODO: emit FIND_THRESHOLD versions


    DEFINE_CLONES(
        ch_repertoire
        .map{ it -> [ it[1] ] }
        .collect(),
        clone_threshold,
        ch_imgt
    )

    emit:
    DEFINE_CLONES.out.tab
}
