include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT } from '../../modules/local/unzip_db'
include { CHANGEO_ASSIGNGENES } from '../../modules/local/changeo/changeo_assigngenes'

workflow VDJ_ANNOTATION {

    take:
    ch_fasta

    main:
    ch_versions = Channel.empty()

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

    CHANGEO_ASSIGNGENES (
        ch_fasta,
        ch_igblast.collect()
    )
    // TODO: check what this does
    //ch_file_sizes = ch_file_sizes.mix(CHANGEO_ASSIGNGENES_REVEAL.out.logs)
    ch_versions = ch_versions.mix(CHANGEO_ASSIGNGENES.out.versions.ifEmpty(null))

    emit:
    fasta_assigned = CHANGEO_ASSIGNGENES.out.fasta
    versions = ch_versions

}
