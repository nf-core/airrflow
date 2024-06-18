include { MIXCR_IND_POSTANALYSIS        } from '../../modules/local/mixcr/mixcr_individualpostanalysis'
include { MIXCR_IND_PLOTS               } from '../../modules/local/mixcr/mixcr_individualpostanalysis_plots'
include { MIXCR_OVERLAP_POSTANALYSIS    } from '../../modules/local/mixcr/mixcr_overlappostanalysis'
include { MIXCR_OVERLAP_PLOTS           } from '../../modules/local/mixcr/mixcr_overlappostanalysis_plots'


workflow MIXCR_POSTANALYSIS {

    take:
    ch_mixcr_clns

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    ch_mixcr_clns.map{ it -> [ it[0].subject_id,
                            it[0].id,
                            it[0].species,
                            it[0].single_cell,
                            it[0].locus,
                            it[1] ] }
            .groupTuple()
            .map{ get_meta_tabs(it) }
            .set { ch_clns_per_subject }
    
            

    MIXCR_IND_POSTANALYSIS (
        ch_clns_per_subject,
        params.mixcr_downsampling,
        params.mixcr_weightfunction,
        params.mixcr_productive_only,
        params.mixcr_drop_outliers
    )

    ch_mixcr_ind_json = MIXCR_IND_POSTANALYSIS.out.mixcr_ind_json

    MIXCR_IND_PLOTS (
        ch_mixcr_ind_json,
        params.mixcr_diversity_plottype,
        params.mixcr_statistical_method,
        params.mixcr_p_adjust_method
    )

    MIXCR_OVERLAP_POSTANALYSIS (
        ch_clns_per_subject,
        params.mixcr_downsampling,
        params.mixcr_weightfunction,
        params.mixcr_productive_only,
        params.mixcr_drop_outliers,
        params.mixcr_overlap_criteria
    )

    ch_mixcr_overlap_json = MIXCR_OVERLAP_POSTANALYSIS.out.mixcr_overlap_json

    ch_mixcr_overlap_json.view()

    MIXCR_OVERLAP_PLOTS {
        ch_mixcr_overlap_json
    }



}


// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id            = [arr[0]].unique().join("")
    meta.sample_ids         = arr[1]
    meta.species            = arr[2]
    meta.single_cell        = arr[3].unique().join("")
    meta.locus              = arr[4].unique().join("")

    def array = []

        array = [ meta, arr[5].flatten() ]

    return array
}