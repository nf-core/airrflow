/*
 * Get tables and group together the ones from the same subject
 */

params.options = [:]

include { MERGE_TABLES } from '../../modules/local/merge_tables' addParams( options: params.options )

workflow MERGE_TABLES_WF {
    take:
    tables

    main:
    tables
        .dump()
        .map{it -> [ it[0].subject+'_'+it[0].locus, it[0].id, it[0].locus, it[0].subject, it[0].species, it[1] ]}
        .groupTuple()
        .dump()
        .map{ get_meta_tabs(it) }
        .dump()
        .set{ch_merge_tables}

    MERGE_TABLES( ch_merge_tables )

    emit:
    tab = MERGE_TABLES.out.tab // channel: [ val(meta), tab ]
}

// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id           = arr[0]
    meta.samples      = arr[1]
    meta.locus        = arr[2].unique().join("")
    meta.subject      = arr[3].unique().join("")
    meta.species      = arr[4].unique().join("")

    def array = []

        array = [ meta, arr[5].flatten() ]

    return array
}
