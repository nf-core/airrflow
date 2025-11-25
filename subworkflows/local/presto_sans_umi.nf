// Include statements

include { GUNZIP            as GUNZIP_SANS_UMI                } from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_SANS_UMI } from '../../modules/local/fastqc_postassembly'
include { FASTP                                               } from '../../modules/nf-core/fastp/main'

//PRESTO
include { PRESTO_ASSEMBLEPAIRS               as  PRESTO_ASSEMBLEPAIRS_SANS_UMI }               from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_FILTERSEQ_POSTASSEMBLY      as  PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI }      from '../../modules/local/presto/presto_filterseq_postassembly'
include { PRESTO_MASKPRIMERS_SCORE           as  PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD }         from '../../modules/local/presto/presto_maskprimers_score'
include { PRESTO_MASKPRIMERS_SCORE           as  PRESTO_MASKPRIMERS_SCORE_SANSUMI_REV }         from '../../modules/local/presto/presto_maskprimers_score'
include { PRESTO_MASKPRIMERS_ALIGN           as  PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD }         from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_ALIGN           as  PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV }         from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_PARSEHEADERS_PRIMERS        as  PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI }         from '../../modules/local/presto/presto_parseheaders_primers'
include { PRESTO_PARSEHEADERS_METADATA       as  PRESTO_PARSEHEADERS_METADATA_SANS_UMI }        from '../../modules/local/presto/presto_parseheaders_metadata'
include { PRESTO_COLLAPSESEQ                 as  PRESTO_COLLAPSESEQ_SANS_UMI }                  from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_SPLITSEQ                    as  PRESTO_SPLITSEQ_SANS_UMI}                      from '../../modules/local/presto/presto_splitseq'


workflow PRESTO_SANS_UMI {
    take:
    ch_reads       // channel: [ val(meta), [ reads ] ]
    ch_cprimers    // channel: [ cprimers.fasta ]
    ch_vprimers    // channel: [ vprimers.fasta ]
    ch_adapter_fasta // channel: [ adapters.fasta ]

    main:

    ch_versions = Channel.empty()

    // Fastp
    save_merged = false
    FASTP (
        ch_reads,
        ch_adapter_fasta,
        params.save_trimmed,
        save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    ch_gunzip = FASTP.out.reads.map{ meta,reads -> [meta, reads[0], reads[1]] }

    // gunzip fastq.gz to fastq
    GUNZIP_SANS_UMI ( ch_gunzip )
    ch_versions = ch_versions.mix(GUNZIP_SANS_UMI.out.versions)

    // Assemble read pairs
    PRESTO_ASSEMBLEPAIRS_SANS_UMI (
        GUNZIP_SANS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_ASSEMBLEPAIRS_SANS_UMI.out.versions)

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI (
        PRESTO_ASSEMBLEPAIRS_SANS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.versions)

    // Mask primers
    def suffix_FWD = "R1"
    def suffix_REV = "R2"
    def barcode_R1 = false
    def barcode_R2 = false
    ch_reads = PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.reads
    if (params.maskprimers_align){
        if (params.cprimer_position == "R1") {
            PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD(
                ch_reads,
                ch_cprimers.collect(),
                params.primer_maxlen,
                params.primer_r1_maxerror,
                params.primer_r1_mask_mode,
                suffix_FWD
            )
            PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV(
                PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD.out.reads,
                ch_vprimers.collect(),
                params.primer_maxlen,
                params.primer_r2_maxerror,
                params.primer_r2_mask_mode,
                suffix_REV
            )
        } else if (params.cprimer_position == "R2") {
            PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD(
                ch_reads,
                ch_vprimers.collect(),
                params.primer_maxlen,
                params.primer_r1_maxerror,
                params.primer_r1_mask_mode,
                suffix_FWD
            )
            PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV(
                PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD.out.reads,
                ch_cprimers.collect(),
                params.primer_maxlen,
                params.primer_r2_maxerror,
                params.primer_r2_mask_mode,
                suffix_REV
            )
        } else {
            error "Error in determining cprimer position. Please choose R1 or R2."
        }

        ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD.out.versions)
        ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV.out.versions)

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_ALIGN_SANSUMI_FWD.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV.out.logs)

        ch_masked_reads = PRESTO_MASKPRIMERS_ALIGN_SANSUMI_REV.out.reads
    } else {
        if (params.cprimer_position == "R1") {
            def start_FWD = "${params.cprimer_start}"
            def start_REV = "${params.vprimer_start}"
            def revpr_FWD = false
            PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD(
                ch_reads,
                ch_cprimers.collect(),
                start_FWD,
                barcode_R1,
                params.primer_r1_maxerror,
                params.primer_r1_mask_mode,
                revpr_FWD,
                suffix_FWD
            )
            PRESTO_MASKPRIMERS_SCORE_SANSUMI_REV(
                PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD.out.reads,
                ch_vprimers.collect(),
                start_REV,
                barcode_R2,
                params.primer_r2_maxerror,
                params.primer_r2_mask_mode,
                params.primer_revpr,
                suffix_REV
            )
        } else if (params.cprimer_position == "R2") {
            def start_FWD = "${params.vprimer_start}"
            def start_REV = "${params.cprimer_start}"
            def revpr_FWD = false

            PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD(
                ch_reads,
                ch_vprimers.collect(),
                start_FWD,
                barcode_R1,
                params.primer_r1_maxerror,
                params.primer_r1_mask_mode,
                revpr_FWD,
                suffix_FWD
            )
            PRESTO_MASKPRIMERS_SCORE_SANSUMI_REV(
                PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD.out.reads,
                ch_cprimers.collect(),
                start_REV,
                barcode_R2,
                params.primer_r2_maxerror,
                params.primer_r2_mask_mode,
                params.primer_revpr,
                suffix_REV
            )
        } else {
            error "Error in determining cprimer position. Please choose R1 or R2."
        }

        ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD.out.versions)

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_SCORE_SANSUMI_FWD.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_SCORE_SANSUMI_REV.out.logs)

        ch_masked_reads = PRESTO_MASKPRIMERS_SCORE_SANSUMI_REV.out.reads

    }

    // Generate QC stats after reads paired and filtered but before collapsed
    FASTQC_POSTASSEMBLY_SANS_UMI (
        ch_masked_reads
    )
    ch_versions = ch_versions.mix(FASTQC_POSTASSEMBLY_SANS_UMI.out.versions)

    // Annotate primers in C_PRIMER and V_PRIMER field
    PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI (
        ch_masked_reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI.out.versions)

    // Annotate metadata on primer headers
    PRESTO_PARSEHEADERS_METADATA_SANS_UMI (
        PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_METADATA_SANS_UMI.out.versions)

    // Mark and count duplicate sequences (DUPCOUNT)
    PRESTO_COLLAPSESEQ_SANS_UMI (
        PRESTO_PARSEHEADERS_METADATA_SANS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_COLLAPSESEQ_SANS_UMI.out.versions)

    // Filter out sequences with less than 2 representative duplicates
    PRESTO_SPLITSEQ_SANS_UMI (
        PRESTO_COLLAPSESEQ_SANS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_SPLITSEQ_SANS_UMI.out.versions)

    emit:
    fasta = PRESTO_SPLITSEQ_SANS_UMI.out.fasta
    versions = ch_versions
    fastp_reads_json = FASTP.out.json.collect{ meta,json -> json }
    fastp_reads_html = FASTP.out.html.collect{ meta,html -> html }
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_SANS_UMI.out.zip
    presto_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_SANS_UMI.out.logs.collect()
    presto_filterseq_logs = PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.logs
    presto_maskprimers_logs = ch_maskprimers_logs.collect()
    presto_collapseseq_logs = PRESTO_COLLAPSESEQ_SANS_UMI.out.logs.collect()
    presto_splitseq_logs = PRESTO_SPLITSEQ_SANS_UMI.out.logs.collect()
}
