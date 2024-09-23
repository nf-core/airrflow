include { FETCH_DATABASES                                               } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST                                     } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT                                        } from '../../modules/local/unzip_db'
include { MIXCR_MIXCR                                                   } from '../../modules/local/mixcr/mixcr'
include { MIXCR_MIXCREXPORTAIRR                                         } from '../../modules/local/mixcr/mixcr_exportairr'
include { MIXCR_MIXCRQCALIGN                                            } from '../../modules/local/mixcr/mixcr_qc_align'
include { MIXCR_MIXCRQCCOVERAGE                                         } from '../../modules/local/mixcr/mixcr_qc_coverage'
include { MIXCR_MIXCRQCTAGS                                             } from '../../modules/local/mixcr/mixcr_qc_tags'
include { MIXCR_MIXCRQCCHAINUSAGE                                       } from '../../modules/local/mixcr/mixcr_qc_chainusage'
include { FASTQ_INPUT_CHECK                                             } from '../../subworkflows/local/fastq_input_check'
include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR  } from '../../modules/local/changeo/changeo_convertdb_fasta'


workflow MIXCR_FLOW {

    take:
    ch_input

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    //
    // read in samplesheet, validate and stage input fies
    //
    FASTQ_INPUT_CHECK(
        ch_input
    )
    ch_versions = ch_versions.mix(FASTQ_INPUT_CHECK.out.versions)

    ch_reads = FASTQ_INPUT_CHECK.out.reads

    MIXCR_MIXCR (
        ch_reads,
        file(params.imgt_mixcr),
        params.kit
    )
    ch_versions = ch_versions.mix(MIXCR_MIXCR.out.versions.first())

    ch_mixcr_out = MIXCR_MIXCR.out.outs


    MIXCR_MIXCREXPORTAIRR (
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr) // it doesnt directly use the imgt db, but it needs it in the right directory anyway
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCREXPORTAIRR.out.versions.first())

    ch_mixcr_airr = MIXCR_MIXCREXPORTAIRR.out.mixcr_airr

    // QC
    MIXCR_MIXCRQCALIGN (
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr)
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCALIGN.out.versions.first())

    MIXCR_MIXCRQCCOVERAGE (
        MIXCR_MIXCR.out.vdjca,
        file(params.imgt_mixcr)
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCCOVERAGE.out.versions.first())

    MIXCR_MIXCRQCCHAINUSAGE (
        MIXCR_MIXCR.out.clns,
        file(params.imgt_mixcr)
        )
    ch_versions = ch_versions.mix(MIXCR_MIXCRQCCHAINUSAGE.out.versions.first())


    // convert airr tsv to fasta
    CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                MIXCR_MIXCREXPORTAIRR.out.mixcr_airr
            )

    ch_versions = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.versions

    ch_fasta = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta

    emit:
    versions = ch_versions
    outs = ch_mixcr_out
    // mixcr output in airr format
    airr = ch_mixcr_airr
    // mixcr output in clns format
    clns = MIXCR_MIXCR.out.clns
    // mixcr output converted to FASTA format
    fasta = ch_fasta
    samplesheet = FASTQ_INPUT_CHECK.out.samplesheet

}

