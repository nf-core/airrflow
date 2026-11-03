include { FIND_THRESHOLD as FIND_CLONAL_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { FIND_THRESHOLD as REPORT_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { CLONAL_ASSIGNMENT } from '../../modules/local/enchantr/clonal_assignment'
include { REPERTOIRE_ANALYSIS} from '../../modules/local/enchantr/repertoire_analysis'
include { DOWSER_LINEAGES } from '../../modules/local/enchantr/dowser_lineages'

workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire_reference
    ch_logo
    clonal_threshold
    skip_report_threshold
    cloneby
    skip_all_clones_report
    lineage_trees
    genotypeby
    crossby
    singlecell
    lineage_tree_builder
    lineage_tree_exec

    main:
    ch_logs = channel.empty()

    if (clonal_threshold == 'auto') {

        ch_find_threshold = ch_repertoire_reference.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        FIND_CLONAL_THRESHOLD (
            ch_find_threshold,
            ch_logo,
            ch_find_threshold_samplesheet,
            cloneby,
            crossby,
            singlecell
        )
        def ch_threshold = FIND_CLONAL_THRESHOLD.out.mean_threshold

        // Collect raw threshold values into a single list so we can distinguish
        // between (A) no values at all (likely upstream failure), and
        // (B) values present but all invalid thresholds ('' / 'NA' / 'NaN').
        def raw_list = ch_threshold
            .splitText( limit:1 ) { it.trim().toString() }
            .map { it -> it.trim() }
            .collect()

        // Process the collected list to identify when no valid thresholds were found
        clone_threshold = raw_list
            .map { list ->
                if (!list || list.size() == 0) {
                    // upstream produced nothing — do not print a message here
                    return []
                }

                def valid = list.findAll { it != '' && it != 'NA' && it != 'NaN' }
                if (valid.size() == 0) {
                    // The automatic threshold finder returned values but all were
                    // NA, NaN or empty strings - ask the user to set a manual value.
                    error "Automatic clone_threshold detection failed. Consider setting --clonal_threshold manually."
                }

                return valid
            }
            .flatten()

    } else {
        clone_threshold = clonal_threshold

        ch_find_threshold = ch_repertoire_reference.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        if (!skip_report_threshold){
            REPORT_THRESHOLD (
                ch_find_threshold,
                ch_logo,
                ch_find_threshold_samplesheet,
                cloneby,
                crossby,
                singlecell
            )
        }
    }

    // Group by the cloneby field and the effective locus: two loci of one subject
    // carry different references and cannot share a clone group.
    ch_repertoire_reference.map{ it ->
                            def meta = it[0]
                            def grouping_locus = meta.locus.toUpperCase()
                            [ meta[cloneby],
                                grouping_locus,
                                meta.id,
                                meta.sample_id,
                                meta.subject_id,
                                meta.species,
                                meta.single_cell,
                                meta.locus,
                                it[1],
                                it[2] ] }
                .set{ ch_repertoire_flat }

    // Only a cloneby value spanning more than one locus needs the suffix, otherwise
    // every existing clone-group output would be renamed.
    ch_repertoire_flat
                .map{ it -> [ it[0], it[1] ] }
                .groupTuple()
                .map{ clone_id, loci -> [ clone_id, loci.unique().size() > 1 ] }
                .set{ ch_mixed_locus }

    ch_repertoire_flat
                .map{ it -> [ it[0], it ] }
                .combine( ch_mixed_locus, by: 0 )
                .map{ _clone_id, row, mixed -> row + [ mixed ] }
                .groupTuple(by: [0,1])
                .map{ get_meta_tabs(it, genotypeby, cloneby) }
                .set{ ch_repertoire_grouped }

    ch_repertoire_grouped.dump(tag: "ch_repertoire_grouped")

    CLONAL_ASSIGNMENT(
        ch_repertoire_grouped,
        clone_threshold.collect(),
        [],
        cloneby,
        singlecell
    )


    // prepare ch for define clones all samples report
    CLONAL_ASSIGNMENT.out.tab
            .map { it -> it[1]}
            .collect()
            .map { it -> [ [id:'all_reps'], it ] }
            .set{ch_all_repertoires_cloned}

    if (!skip_all_clones_report){

        ch_all_repertoires_cloned_samplesheet = ch_all_repertoires_cloned.map{ it -> it[1] }
                                        .collect()
                                        .flatten()
                                        .map{ it -> it.getName().toString() }
                                        .collectFile(name: 'all_repertoires_cloned_samplesheet.txt', newLine: true)

        REPERTOIRE_ANALYSIS(
            ch_all_repertoires_cloned,
            ch_all_repertoires_cloned_samplesheet,
            cloneby
        )
    }

    if (lineage_trees){
        DOWSER_LINEAGES(
            CLONAL_ASSIGNMENT.out.tab,
            lineage_tree_builder,
            lineage_tree_exec
        )
    }

    emit:
    repertoire = REPERTOIRE_ANALYSIS.out.tab
    logs = ch_logs
}

// Function to map
// arr[0] cloneby value, arr[1] grouping locus, arr[10] whether it spans >1 locus.
def get_meta_tabs(arr, genotypeby, cloneby) {
    def clone_id = arr[0]
    def grouping_locus = arr[1]

    if (arr[4].unique().size() > 1) {
            error "Multiple subject_id found for ${clone_id} (${arr[4].join(', ')}). Please check your input parameters and ensure that all samples with the same 'cloneby' value have the same 'subject_id' value."
    }

    def locus_class = arr[7].unique().join("")

    def meta = [:]
    meta.id                 = arr[10].unique().contains(true) ? "${clone_id}_${grouping_locus}".toString() : clone_id.toString()
    meta.sample_id          = arr[3].flatten()
    meta.subject_id         = arr[4].unique().join("")
    meta.species            = arr[5].unique().join("")
    meta.single_cell        = arr[6].unique().join("")
    meta.locus              = locus_class

    def array = []

        array = [ meta, arr[8].flatten(), arr[9].unique() ]
        if (arr[9].size() > 1) {
            error "Multiple reference fasta files found for ${meta.id}. Please check your input parameters and ensure that all samples with the same ${genotypeby} value (parameter 'genotype_by') have the same ${cloneby} value (parameter 'clone_by')."
        }
    return array
}
