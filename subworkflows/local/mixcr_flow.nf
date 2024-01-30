// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FETCH_DATABASES               } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST     } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT        } from '../../modules/local/unzip_db'
include { MIXCR_MIXCR                   } from '../../modules/local/mixcr/mixcr'
include { MIXCR_MIXCREXPORTAIRR         } from '../../modules/local/mixcr/mixcr_exportairr'
include { MIXCR_MIXCRQCALIGN            } from '../../modules/local/mixcr/mixcr_qc_align'
include { MIXCR_MIXCRQCCOVERAGE         } from '../../modules/local/mixcr/mixcr_qc_coverage'
include { MIXCR_MIXCRQCTAGS             } from '../../modules/local/mixcr/mixcr_qc_tags'
include { MIXCR_MIXCRQCCHAINUSAGE       } from '../../modules/local/mixcr/mixcr_qc_chainusage'

include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR  } from '../../modules/local/changeo/changeo_convertdb_fasta'


workflow MIXCR_FLOW {

    take:
    ch_reads_bulk // meta, reads

    main:

    ch_reads_bulk.dump(tag: "ch_reads_bulk")

    ch_versions = Channel.empty()

    // // FETCH DATABASES
    // // TODO: this can take a long time, and the progress shows 0%. Would be
    // // nice to have some better progress reporting.
    // // And maybe run this as 2 separate steps, one for IMGT and one for IgBLAST?
    // if( !params.fetch_imgt ){
    //     if (params.igblast_base.endsWith(".zip")) {
    //         Channel.fromPath("${params.igblast_base}")
    //                 .ifEmpty{ error "IGBLAST DB not found: ${params.igblast_base}" }
    //                 .set { ch_igblast_zipped }
    //         UNZIP_IGBLAST( ch_igblast_zipped.collect() )
    //         ch_igblast = UNZIP_IGBLAST.out.unzipped
    //         ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions.ifEmpty(null))
    //     } else {
    //         Channel.fromPath("${params.igblast_base}")
    //             .ifEmpty { error "IGBLAST DB not found: ${params.igblast_base}" }
    //             .set { ch_igblast }
    //     }
    // }

    // if( !params.fetch_imgt ){
    //     if (params.imgtdb_base.endsWith(".zip")) {
    //         Channel.fromPath("${params.imgtdb_base}")
    //                 .ifEmpty{ error "IMGTDB not found: ${params.imgtdb_base}" }
    //                 .set { ch_imgt_zipped }
    //         UNZIP_IMGT( ch_imgt_zipped.collect() )
    //         ch_imgt = UNZIP_IMGT.out.unzipped
    //         ch_versions = ch_versions.mix(UNZIP_IMGT.out.versions.ifEmpty(null))
    //     } else {
    //         Channel.fromPath("${params.imgtdb_base}")
    //             .ifEmpty { error "IMGT DB not found: ${params.imgtdb_base}" }
    //             .set { ch_imgt }
    //     }
    // }

    // if (params.fetch_imgt) {
    //     FETCH_DATABASES()
    //     ch_igblast = FETCH_DATABASES.out.igblast
    //     ch_imgt = FETCH_DATABASES.out.imgt
    //     ch_versions = ch_versions.mix(FETCH_DATABASES.out.versions.ifEmpty(null))
    // }

    MIXCR_MIXCR ( 
        ch_reads_bulk,
        file(params.imgt_mixcr),
        params.kit
    )
    ch_versions = ch_versions.mix(MIXCR_MIXCR.out.versions.first())


    MIXCR_MIXCREXPORTAIRR ( 
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr) // it doesnt directly use the imgt db, but it needs it in the right directory anyway
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCREXPORTAIRR.out.versions.first())

    // QC
    MIXCR_MIXCRQCALIGN ( 
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr) // it doesnt directly use the imgt db, but it needs it in the right directory anyway
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCALIGN.out.versions.first())

    MIXCR_MIXCRQCCOVERAGE ( 
        MIXCR_MIXCR.out.vdjca,
        file(params.imgt_mixcr) // it doesnt directly use the imgt db, but it needs it in the right directory anyway
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCCOVERAGE.out.versions.first())

    MIXCR_MIXCRQCCHAINUSAGE ( 
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr) // it doesnt directly use the imgt db, but it needs it in the right directory anyway
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCCHAINUSAGE.out.versions.first())


    // convert airr tsv to fasta (cellranger does not create any fasta with clonotype information)
    CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                MIXCR_MIXCREXPORTAIRR.out.mixcr_airr
            )

    
    ch_fasta = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta

    emit:
    versions = ch_versions
    // complete cellranger output
    fasta = ch_fasta
}

