/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { MERGE_TABLES } from '../modules/local/merge_tables' addParams( options: params.options )

workflow MERGE_TABLES_WF {
    take:
    tables // file: /path/to/samplesheet.tsv
    
    main:
    // TODO: avoiding checking samplesheet for now, add samplesheet check later.
    tables        
        .dump()
        .map{it -> [ it[0].source, it[0].id, it[1] ]}
        .groupTuple()
        .dump()
        .map{ get_meta_tabs(it) }
        .dump()
        .set{ch_merge_tables}
    
    MERGE_TABLES( ch_merge_tables )

    emit:
    MERGE_TABLES.out.tab // channel: [ val(meta), tab ]
}

// Function to map 
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id           = arr[0]
    meta.samples      = arr[1]

    def array = []

        array = [ meta, arr[2].flatten() ]

    return array    
}