include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'

workflow DATABASES {

    take:
    fetch_germlines
    reference_igblast
    reference_fasta

    main:

    if( fetch_germlines == "none"){
        if (reference_igblast.endsWith(".zip")) {
            channel.fromPath("${reference_igblast}")
                    .ifEmpty{ error "IGBLAST DB not found: ${reference_igblast}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
        } else {
            channel.fromPath("${reference_igblast}")
                .ifEmpty { error "IGBLAST DB not found: ${reference_igblast}" }
                .set { ch_igblast }
        }
    }

    if( fetch_germlines == "none" ){
        if (reference_fasta.endsWith(".zip")) {
            channel.fromPath("${reference_fasta}")
                    .ifEmpty{ error "IMGTDB not found: ${reference_fasta}" }
                    .set { ch_reference_fasta_zipped }
            UNZIP_REFERENCE_FASTA( ch_reference_fasta_zipped.collect() )
            ch_reference_fasta = UNZIP_REFERENCE_FASTA.out.unzipped
        } else {
            channel.fromPath("${reference_fasta}")
                .ifEmpty { error "IMGT DB not found: ${reference_fasta}" }
                .set { ch_reference_fasta }
        }

        VALIDATE_IGBLAST_DB(ch_igblast, ch_reference_fasta)
        ch_igblast = VALIDATE_IGBLAST_DB.out.igblast
    }

    if (fetch_germlines == "imgt" || fetch_germlines == "airrc-imgt") {
        FETCH_DATABASES(channel.value(fetch_germlines))
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_reference_fasta = FETCH_DATABASES.out.reference_fasta
    }

    emit:
    reference_fasta = ch_reference_fasta
    igblast = ch_igblast
}
