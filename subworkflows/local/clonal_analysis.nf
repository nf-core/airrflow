include { FIND_THRESHOLD as FIND_CLONAL_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { FIND_THRESHOLD as REPORT_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { CLONAL_ASSIGNMENT as CLONAL_ASSIGNMENT_COMPUTE  } from '../../modules/local/enchantr/clonal_assignment'
include { CLONAL_ASSIGNMENT as CLONAL_ASSIGNMENT_REPORT } from '../../modules/local/enchantr/clonal_assignment'
include { DOWSER_LINEAGES } from '../../modules/local/enchantr/dowser_lineages'

workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire
    ch_reference_fasta
    ch_logo

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    if (params.clonal_threshold == "auto") {

        ch_find_threshold = ch_repertoire.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .dump(tag: 'ch_find_threshold_samplesheet')
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        FIND_CLONAL_THRESHOLD (
            ch_find_threshold,
            ch_logo,
            ch_find_threshold_samplesheet
        )
        ch_threshold = FIND_CLONAL_THRESHOLD.out.mean_threshold
        ch_versions = ch_versions.mix(FIND_CLONAL_THRESHOLD.out.versions)

        clone_threshold = ch_threshold
            .splitText( limit:1 ) { it.trim().toString() }
            .dump(tag: 'clone_threshold')
            .filter { it != 'NA'}
            .filter { it != 'NaN' }
            .ifEmpty { error "Automatic clone_threshold is 'NA'. Consider setting --clonal_threshold manually."}

    } else {
        clone_threshold = params.clonal_threshold

        ch_find_threshold = ch_repertoire.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .dump(tag: 'ch_find_threshold_samplesheet')
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        if (!params.skip_report_threshold){
            REPORT_THRESHOLD (
                ch_find_threshold,
                ch_logo,
                ch_find_threshold_samplesheet
            )
            ch_versions = ch_versions.mix(REPORT_THRESHOLD.out.versions)
        }
    }

    // prepare ch for define clones
    ch_repertoire.map{ it -> [ it[0]."${params.cloneby}",
                                it[0].id,
                                it[0].subject_id,
                                it[0].species,
                                it[0].single_cell,
                                it[0].locus,
                                it[1] ] }
                .groupTuple()
                .map{ get_meta_tabs(it) }
                .set{ ch_clonal_assignment }

    CLONAL_ASSIGNMENT_COMPUTE(
        ch_clonal_assignment,
        clone_threshold.collect(),
        ch_reference_fasta.collect(),
        []
    )

    ch_versions = ch_versions.mix(CLONAL_ASSIGNMENT_COMPUTE.out.versions)

    // prepare ch for define clones all samples report
    CLONAL_ASSIGNMENT_COMPUTE.out.tab
            .map { it -> it[1]}
            .collect()
            .map { it -> [ [id:'all_reps'], it ] }
            .set{ch_all_repertoires_cloned}

    if (!params.skip_all_clones_report){

        ch_all_repertoires_cloned_samplesheet = ch_all_repertoires_cloned.map{ it -> it[1] }
                                        .collect()
                                        .flatten()
                                        .map{ it -> it.getName().toString() }
                                        .dump(tag: 'ch_all_repertoires_cloned_samplesheet')
                                        .collectFile(name: 'all_repertoires_cloned_samplesheet.txt', newLine: true)

        CLONAL_ASSIGNMENT_REPORT(
            ch_all_repertoires_cloned,
            clone_threshold.collect(),
            ch_reference_fasta.collect(),
            ch_all_repertoires_cloned_samplesheet
        )
        ch_versions = ch_versions.mix(CLONAL_ASSIGNMENT_REPORT.out.versions)
    }

    if (params.lineage_trees){
        DOWSER_LINEAGES(
            CLONAL_ASSIGNMENT_COMPUTE.out.tab
        )
        ch_versions = ch_versions.mix(DOWSER_LINEAGES.out.versions)
    }

    emit:
    repertoire = CLONAL_ASSIGNMENT_COMPUTE.out.tab
    versions = ch_versions
    logs = ch_logs
}

// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id            = [arr[0]].unique().join("")
    meta.sample_ids         = arr[1]
    meta.subject_id         = arr[2]
    meta.species            = arr[3]
    meta.single_cell        = arr[4].unique().join("")
    meta.locus              = arr[5].unique().join("")

    def array = []

        array = [ meta, arr[6].flatten() ]

    return array
}
