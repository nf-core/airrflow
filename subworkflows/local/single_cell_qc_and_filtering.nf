include { SINGLE_CELL_QC  } from '../../modules/local/enchantr/single_cell_qc'

workflow SINGLE_CELL_QC_AND_FILTERING {
    take:
    repertoires // tuple [meta, repertoire_tab]

    main:
    ch_versions = Channel.empty()

    repertoires
            .dump(tag:"scqc-reps")
            .map{ it -> [   it[0].single_cell,
                            it[0].id,
                            it[0].filename,
                            it[0].subject_id,
                            it[0].species,
                            it[0].collapseby_group,
                            it[0].collapseby_size,
                            it[0].cloneby_group,
                            it[0].cloneby_size,
                            it[0].filetype,
                            it[0].pcr_target_locus,
                            it[0].locus,
                            it[1]  ] }
            .groupTuple()
            .dump(tag:'scqc-aftergroup')
            .map{ get_meta_tabs(it) }
            .dump(tag: 'scqc-aftergetmeta')
            .set{ch_repertoire_all}

    repertoires
        .collectFile( keepHeader: true, skip:1, storeDir: "${params.outdir}/csv") { meta, tab ->
                        sample_id   = meta.id
                        orig_filename    = meta.filename
                        subject_id  = meta.subject_id
                        species     = meta.species
                        collapseby_group    = meta.collapseby_group
                        collapseby_size     = meta.collapseby_size
                        cloneby_group       = meta.cloneby_group
                        cloneby_size        = meta.cloneby_size
                        filetype            = meta.filetype
                        pcr_target_locus    = meta.pcr_target_locus
                        locus               = meta.locus
                        tab                 = "${params.outdir}/qc-filtering/single-cell-qc/all_reps_scqc_report/${sample_id}__scqc-pass.tsv"

                        ["single_cell_qc.csv", "sample_id,orig_filename,subject_id,species,collapseby_group,collapseby_size,cloneby_group,cloneby_size,filetype,pcr_target_locus,locus,tab\n${sample_id},${orig_filename},${subject_id},${species},${collapseby_group},${collapseby_size},${cloneby_group},${cloneby_size},${filetype},${pcr_target_locus},${locus},${tab}\n"]
                    }

    SINGLE_CELL_QC(
        ch_repertoire_all
    )
    // ch_file_sizes = ch_file_sizes.mix(SINGLE_CELL_QC.out.logs)
    ch_versions = ch_versions.mix(SINGLE_CELL_QC.out.versions.ifEmpty(null))

    emit:
    versions = ch_versions
    repertoires = SINGLE_CELL_QC.out.tab
}

// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id                     = arr[1]
    meta.filename               = arr[2]
    meta.subject_id             = arr[3]
    meta.species                = arr[4]
    meta.collapseby_group       = arr[5]
    meta.collapseby_size        = arr[6]
    meta.cloneby_group          = arr[7]
    meta.cloneby_size           = arr[8]
    meta.filetype               = arr[9]
    meta.pcr_target_locus       = arr[10]
    meta.locus                  = arr[11]

    def array = []
        array = [ meta, arr[12].flatten() ]
    return array
}