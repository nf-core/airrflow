include { FETCH_IMGT_REFERENCE } from '../../modules/local/fetch_imgt_reference'
include { FETCH_AIRRC_IMGT_REFERENCE } from '../../modules/local/fetch_airrc_imgt_reference'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'

workflow DATABASES {

    take:

    main:
    ch_versions = Channel.empty()

    // FETCH DATABASES
    if( !params.fetch_reference ){
        if (params.reference_igblast.endsWith(".zip")) {
            Channel.fromPath("${params.reference_igblast}")
                    .ifEmpty{ error "IGBLAST DB not found: ${params.reference_igblast}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions)
        } else {
            Channel.fromPath("${params.reference_igblast}")
                .ifEmpty { error "IGBLAST DB not found: ${params.reference_igblast}" }
                .set { ch_igblast }
        }
    }

    if( !params.fetch_reference ){
        if (params.reference_fasta.endsWith(".zip")) {
            Channel.fromPath("${params.reference_fasta}")
                    .ifEmpty{ error "IMGTDB not found: ${params.reference_fasta}" }
                    .set { ch_reference_fasta_zipped }
            UNZIP_REFERENCE_FASTA( ch_reference_fasta_zipped.collect() )
            ch_reference_fasta = UNZIP_REFERENCE_FASTA.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_REFERENCE_FASTA.out.versions)
        } else {
            Channel.fromPath("${params.reference_fasta}")
                .ifEmpty { error "IMGT DB not found: ${params.reference_fasta}" }
                .set { ch_reference_fasta }
        }
    }

    if (params.fetch_reference == 'imgt') {
        FETCH_IMGT_REFERENCE()
        ch_igblast = FETCH_IMGT_REFERENCE.out.igblast
        ch_reference_fasta = FETCH_IMGT_REFERENCE.out.reference_fasta
        ch_versions = ch_versions.mix(FETCH_IMGT_REFERENCE.out.versions)
    } else if (params.fetch_reference == 'airrc-imgt') {
        FETCH_AIRRC_IMGT_REFERENCE()
        ch_igblast = FETCH_AIRRC_IMGT_REFERENCE.out.igblast
        ch_reference_fasta = FETCH_AIRRC_IMGT_REFERENCE.out.reference_fasta
        ch_versions = ch_versions.mix(FETCH_AIRRC_IMGT_REFERENCE.out.versions)
    }

    emit:
    versions = ch_versions
    reference_fasta = ch_reference_fasta
    igblast = ch_igblast
}
