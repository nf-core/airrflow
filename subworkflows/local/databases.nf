include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT } from '../../modules/local/unzip_db'

workflow DATABASES {

    take:

    main:
    ch_versions = Channel.empty()

    // FETCH DATABASES
    if( !params.fetch_imgt ){
        if (params.igblast_base.endsWith(".zip")) {
            Channel.fromPath("${params.igblast_base}")
                    .ifEmpty{ error "IGBLAST DB not found: ${params.igblast_base}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions)
        } else {
            Channel.fromPath("${params.igblast_base}")
                .ifEmpty { error "IGBLAST DB not found: ${params.igblast_base}" }
                .set { ch_igblast }
        }
    }

    if( !params.fetch_imgt ){
        if (params.imgtdb_base.endsWith(".zip")) {
            Channel.fromPath("${params.imgtdb_base}")
                    .ifEmpty{ error "IMGTDB not found: ${params.imgtdb_base}" }
                    .set { ch_imgt_zipped }
            UNZIP_IMGT( ch_imgt_zipped.collect() )
            ch_imgt = UNZIP_IMGT.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IMGT.out.versions)
        } else {
            Channel.fromPath("${params.imgtdb_base}")
                .ifEmpty { error "IMGT DB not found: ${params.imgtdb_base}" }
                .set { ch_imgt }
        }
    }

    if (params.fetch_imgt) {
        FETCH_DATABASES()
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
        ch_versions = ch_versions.mix(FETCH_DATABASES.out.versions)
    }

    emit:
    versions = ch_versions
    imgt = ch_imgt
    igblast = ch_igblast
}
