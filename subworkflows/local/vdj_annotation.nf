include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT } from '../../modules/local/unzip_db'
include { CHANGEO_ASSIGNGENES } from '../../modules/local/changeo/changeo_assigngenes'
include { CHANGEO_MAKEDB } from '../../modules/local/changeo/changeo_makedb'
include { CHANGEO_PARSEDB_SPLIT } from '../../modules/local/changeo/changeo_parsedb_split'
include { IGBLAST_ASSIGNGENES } from '../../modules/local/igblast/igblast_assigngenes'
// reveal
include { FILTER_QUALITY  } from '../../modules/local/reveal/filter_quality'
include { FILTER_JUNCTION_MOD3  } from '../../modules/local/reveal/filter_junction_mod3'
include { ADD_META_TO_TAB  } from '../../modules/local/reveal/add_meta_to_tab'


workflow VDJ_ANNOTATION {

    take:
    ch_fasta // [meta, fasta]
    ch_validated_samplesheet

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // FETCH DATABASES
    // TODO: this can take a long time, and the progress shows 0%. Would be
    // nice to have some better progress reporting.
    // And maybe run this as 2 separate steps, one for IMGT and one for IgBLAST?
    if( params.igblast_base ){
        if (params.igblast_base.endsWith(".zip")) {
            Channel.fromPath("${params.igblast_base}")
                    .ifEmpty{ exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions.ifEmpty(null))
        } else {
            Channel.fromPath("${params.igblast_base}")
                .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
                .set { ch_igblast }
        }
    }

    if( params.imgtdb_base ){
        if (params.imgtdb_base.endsWith(".zip")) {
            Channel.fromPath("${params.imgtdb_base}")
                    .ifEmpty{ exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
                    .set { ch_imgt_zipped }
            UNZIP_IMGT( ch_imgt_zipped.collect() )
            ch_imgt = UNZIP_IMGT.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IMGT.out.versions.ifEmpty(null))
        } else {
            Channel.fromPath("${params.imgtdb_base}")
                .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
                .set { ch_imgt }
        }
    }

    if (!params.igblast_base | !params.imgtdb_base) {
        FETCH_DATABASES()
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
        ch_versions = ch_versions.mix(FETCH_DATABASES.out.versions.ifEmpty(null))
    }
    if (!params.directcall_igblast){
        CHANGEO_ASSIGNGENES (
            ch_fasta,
            ch_igblast.collect()
        )

        ch_logs = ch_logs.mix(CHANGEO_ASSIGNGENES.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_ASSIGNGENES.out.versions.ifEmpty(null))

        CHANGEO_MAKEDB (
            CHANGEO_ASSIGNGENES.out.fasta,
            CHANGEO_ASSIGNGENES.out.blast,
            ch_imgt.collect()
        )
        ch_logs = ch_logs.mix(CHANGEO_MAKEDB.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_MAKEDB.out.versions.ifEmpty(null))

        ch_assigned_tab = CHANGEO_MAKEDB.out.tab
        ch_assignment_logs = CHANGEO_MAKEDB.out.logs

    } else {
        IGBLAST_ASSIGNGENES(
            ch_fasta,
            ch_igblast.collect()
        )

        ch_assigned_tab = IGBLAST_ASSIGNGENES.out.tab

        ch_logs = ch_logs.mix(IGBLAST_ASSIGNGENES.out.logs)
        ch_versions = ch_versions.mix(IGBLAST_ASSIGNGENES.out.versions.ifEmpty(null))
        ch_assignment_logs = IGBLAST_ASSIGNGENES.out.makedb_log
    }

    // Apply quality filters:
    // - locus should match v_call chain
    // - seq alignment min length informative positions 200
    // - max 10% N nucleotides
    // TODO: emit versions
    FILTER_QUALITY(
        ch_assigned_tab
    )
    ch_logs = ch_logs.mix(FILTER_QUALITY.out.logs)

    if (params.productive_only) {
        CHANGEO_PARSEDB_SPLIT (
            FILTER_QUALITY.out.tab
        )
        ch_logs = ch_logs.mix(CHANGEO_PARSEDB_SPLIT.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_PARSEDB_SPLIT.out.versions.ifEmpty(null))

        // Apply filter: junction length multiple of 3
        // TODO: Add to enchantr and emit versions?
        FILTER_JUNCTION_MOD3(
            CHANGEO_PARSEDB_SPLIT.out.tab
        )
        ch_logs = ch_logs.mix(FILTER_JUNCTION_MOD3.out.logs)
        ch_repertoire = FILTER_JUNCTION_MOD3.out.tab.ifEmpty(null)
    } else {
        ch_repertoire = FILTER_QUALITY.out.tab.ifEmpty(null)
    }

    ADD_META_TO_TAB(
        ch_repertoire,
        ch_validated_samplesheet
    )
    //TODO: emit versions
    ch_logs = ch_logs.mix(ADD_META_TO_TAB.out.logs)


    emit:
    versions = ch_versions
    repertoire = ADD_META_TO_TAB.out.tab
    imgt = ch_imgt
    igblast = ch_igblast
    changeo_makedb_logs = ch_assignment_logs
    logs = ch_logs

}
