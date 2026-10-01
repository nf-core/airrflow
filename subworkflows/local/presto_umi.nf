// Include statements
include { MERGE_UMI                                      } from '../../modules/local/merge_UMI'
include { RENAME_FASTQ        as RENAME_FASTQ_UMI        } from '../../modules/local/rename_fastq'
include { GUNZIP              as GUNZIP_UMI              } from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_UMI } from '../../modules/local/fastqc_postassembly'
include { FASTP                                          } from '../../modules/nf-core/fastp/main'


//PRESTO
include { PRESTO_FILTERSEQ      as  PRESTO_FILTERSEQ_UMI              }    from '../../modules/local/presto/presto_filterseq'
include { PRESTO_MASKPRIMERS    as  PRESTO_MASKPRIMERS_UMI            }    from '../../modules/local/presto/presto_maskprimers'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_ALIGN_PRIMERS            }    from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_EXTRACT as PRESTO_MASKPRIMERS_EXTRACT    }    from '../../modules/local/presto/presto_maskprimers_extract'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_ALIGN_CREGION            }    from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_ALIGN_UMI_R1 }    from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_MASKPRIMERS_ALIGN_UMI_R2 }    from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_SCORE as PRESTO_MASKPRIMERS_SCORE_UMI_R1 }    from '../../modules/local/presto/presto_maskprimers_score'
include { PRESTO_MASKPRIMERS_SCORE as PRESTO_MASKPRIMERS_SCORE_UMI_R2 }    from '../../modules/local/presto/presto_maskprimers_score'
include { PRESTO_MASKPRIMERS_EXTRACT as PRESTO_MASKPRIMERS_EXTRACT_R1 }    from '../../modules/local/presto/presto_maskprimers_extract'
include { PRESTO_MASKPRIMERS_EXTRACT as PRESTO_MASKPRIMERS_EXTRACT_R2 }    from '../../modules/local/presto/presto_maskprimers_extract'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_UMI                }    from '../../modules/local/presto/presto_pairseq'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_ALIGN              }    from '../../modules/local/presto/presto_pairseq'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_EXTRACT            }    from '../../modules/local/presto/presto_pairseq'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_CLUSTERSETS        }    from '../../modules/local/presto/presto_pairseq'
include { PRESTO_CLUSTERSETS    as  PRESTO_CLUSTERSETS_UMI            }    from '../../modules/local/presto/presto_clustersets'
include { PRESTO_PARSE_CLUSTER  as  PRESTO_PARSE_CLUSTER_UMI          }    from '../../modules/local/presto/presto_parse_cluster'
include { PRESTO_BUILDCONSENSUS as  PRESTO_BUILDCONSENSUS_UMI         }    from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_BUILDCONSENSUS as PRESTO_BUILDCONSENSUS_ALIGN_RACE   }    from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_BUILDCONSENSUS as PRESTO_BUILDCONSENSUS_EXTRACT      }    from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_POSTCONSENSUS_PAIRSEQ as PRESTO_POSTCONSENSUS_PAIRSEQ_UMI }    from '../../modules/local/presto/presto_postconsensus_pairseq'
include { PRESTO_ASSEMBLEPAIRS  as  PRESTO_ASSEMBLEPAIRS_UMI          }    from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_ASSEMBLEPAIRS_SEQUENTIAL                             }    from '../../modules/local/presto/presto_assemblepairs_sequential'
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_COLLAPSE_UMI  }    from '../../modules/local/presto/presto_parseheaders'
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_CREGION       }    from '../../modules/local/presto/presto_parseheaders'
include { PRESTO_PARSEHEADERS_PRIMERS   as PRESTO_PARSEHEADERS_PRIMERS_UMI  }   from '../../modules/local/presto/presto_parseheaders_primers'
include { PRESTO_PARSEHEADERS_METADATA  as PRESTO_PARSEHEADERS_METADATA_UMI }   from '../../modules/local/presto/presto_parseheaders_metadata'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_UMI             }    from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_ALIGN           }    from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_CREGION         }    from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_SPLITSEQ       as PRESTO_SPLITSEQ_UMI                }    from '../../modules/local/presto/presto_splitseq'


workflow PRESTO_UMI {
    take:
    ch_reads       // channel: [ val(meta), [ reads ] ]
    ch_cprimers    // channel: [ cprimers.fasta ]
    ch_vprimers    // channel: [ vprimers.fasta ]
    ch_adapter_fasta // channel: [ adapters.fasta ]
    ch_internal_cregion // channel: [ internal_cregions.fasta ]
    ch_igblast
    maskprimers_align_race
    umi_position
    cprimer_position
    index_file
    save_trimmed
    primer_maxlen
    primer_r1_maxerror
    primer_r1_mask_mode
    umi_length
    umi_start
    primer_r2_extract_len
    primer_r2_mask_mode
    maskprimers_align
    primer_r2_maxerror
    primer_revpr
    maskprimers_extract
    cprimer_start
    primer_r1_extract_len
    vprimer_start
    cluster_sets
    assemblepairs_sequential
    align_cregion
    cregion_maxlen
    cregion_maxerror
    cregion_mask_mode
    filterseq_q
    buildconsensus_maxerror
    buildconsensus_maxgap
    primer_consensus

    main:
    ch_versions = channel.empty()

    // Validate params
    if (maskprimers_align_race & umi_position == 'R1') {error "The maskprimers align race option is only supported with UMI barcodes in the R2 reads (reads containing V region)."}
    if (maskprimers_align_race & cprimer_position == 'R2') {error "The maskprimers align race option is only supported with Cprimers in the R1 reads (reads containing C region)."}


    // Merge UMI from index file to R1 if provided
    if (index_file) {

        // ch for fastp reads R1 R2
        ch_reads.map{ meta, reads -> [meta, [reads[0], reads[1]]] }
                .set{ ch_reads_R1_R2 }

        // Fastp reads R1 R2
        save_merged = false
        FASTP (
            ch_reads_R1_R2,
            ch_adapter_fasta,
            save_trimmed,
            save_merged
        )
        ch_versions = ch_versions.mix(FASTP.out.versions)
        //ch for merge umi
        ch_meta_R1_R2 = FASTP.out.reads
                                        .map{ meta, reads -> [meta.id, meta, reads[0], reads[1]] }
        ch_meta_index = ch_reads
                                .map{ meta, reads -> [meta.id, meta, reads[2]] }
        ch_meta_R1_R2_index = ch_meta_R1_R2.join( ch_meta_index )
                                            .map{ id, meta1, R1, R2, meta2, index -> [ meta1, R1, R2, index ] }

        MERGE_UMI (
            ch_meta_R1_R2_index,
            umi_start,
            umi_length
        )
        ch_gunzip = MERGE_UMI.out.reads


    } else {

        // Fastp reads
        save_merged = false
        FASTP (
            ch_reads,
            ch_adapter_fasta,
            save_trimmed,
            save_merged
        )
        ch_versions = ch_versions.mix(FASTP.out.versions)
        ch_rename_fastq_umi = FASTP.out.reads.map{ meta,reads -> [meta, reads[0], reads[1]] }

        RENAME_FASTQ_UMI ( ch_rename_fastq_umi )
        ch_gunzip = RENAME_FASTQ_UMI.out.reads

    }

    // gunzip fastq.gz to fastq
    GUNZIP_UMI ( ch_gunzip )

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_UMI (
        GUNZIP_UMI.out.reads,
        filterseq_q
    )

    // Split reads into R1 and R2
    ch_reads_R1 = PRESTO_FILTERSEQ_UMI.out.reads
                                        .map{ reads -> [reads[0], reads[1]] }
    ch_reads_R2 = PRESTO_FILTERSEQ_UMI.out.reads
                                        .map{ reads -> [reads[0], reads[2]] }

    // Mask primers
    if (maskprimers_align_race) {

        PRESTO_ALIGN_PRIMERS(
            ch_reads_R1,
            ch_cprimers.collect(),
            primer_maxlen,
            false,
            primer_r1_maxerror,
            primer_r1_mask_mode,
            false,
            "R1"
        )

        PRESTO_MASKPRIMERS_EXTRACT(
            ch_reads_R2,
            umi_length,
            primer_r2_extract_len,
            primer_r2_mask_mode,
            true,
            "R2"
        )

        // Merge again R1 and R2 by sample ID.
        ch_maskprimers_reads_R1 = PRESTO_ALIGN_PRIMERS.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads_R2 = PRESTO_MASKPRIMERS_EXTRACT.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads = ch_maskprimers_reads_R1.join(ch_maskprimers_reads_R2)
                                                        .map{ it -> [it[1], it[2], it[4]] }

        ch_maskprimers_logs = PRESTO_ALIGN_PRIMERS.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_EXTRACT.out.logs)

        def barcode_position = (umi_position == "R1" | index_file) ? "R1" : "R2"

        PRESTO_PAIRSEQ_ALIGN(
            ch_maskprimers_reads,
            barcode_position
        )
        ch_for_clustersets = PRESTO_PAIRSEQ_ALIGN.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_ALIGN.out.logs

    } else if (maskprimers_align) {

        if (cprimer_position == "R1") {
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false

            PRESTO_MASKPRIMERS_ALIGN_UMI_R1(
                ch_reads_R1,
                ch_cprimers.collect(),
                primer_maxlen,
                barcode_R1,
                primer_r1_maxerror,
                primer_r1_mask_mode,
                false,
                "R1"
            )

            def barcode_R2 = (umi_position == "R2") ? true : false

            PRESTO_MASKPRIMERS_ALIGN_UMI_R2(
                ch_reads_R2,
                ch_vprimers.collect(),
                primer_maxlen,
                barcode_R2,
                primer_r2_maxerror,
                primer_r2_mask_mode,
                primer_revpr,
                "R2"
            )
        } else if (cprimer_position == "R2") {
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false

            PRESTO_MASKPRIMERS_ALIGN_UMI_R1(
                ch_reads_R1,
                ch_vprimers.collect(),
                primer_maxlen,
                barcode_R1,
                primer_r1_maxerror,
                primer_r1_mask_mode,
                false,
                "R1"
            )
            def barcode_R2 = (umi_position == "R2") ? true : false

            PRESTO_MASKPRIMERS_ALIGN_UMI_R2(
                ch_reads_R2,
                ch_cprimers.collect(),
                primer_maxlen,
                barcode_R2,
                primer_r2_maxerror,
                primer_r2_mask_mode,
                primer_revpr,
                "R2"
            )
        } else {
            error "Error in determining cprimer position. Please choose R1 or R2."
        }


        // Merge again R1 and R2 by sample ID.
        ch_maskprimers_reads_R1 = PRESTO_MASKPRIMERS_ALIGN_UMI_R1.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads_R2 = PRESTO_MASKPRIMERS_ALIGN_UMI_R2.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads = ch_maskprimers_reads_R1.join(ch_maskprimers_reads_R2)
                                                        .map{ it -> [it[1], it[2], it[4]] }

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_ALIGN_UMI_R1.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_ALIGN_UMI_R2.out.logs)

        // Pre-consensus pair
        def barcode_position = (umi_position == "R1" | index_file) ? "R1" : "R2"

        PRESTO_PAIRSEQ_UMI (
            ch_maskprimers_reads,
            barcode_position
        )

        ch_for_clustersets = PRESTO_PAIRSEQ_UMI.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_UMI.out.logs

    } else if (maskprimers_extract){

        if (cprimer_position == "R1"){
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false
            def cprimer_start_var = (index_file | umi_position == 'R1') ? "${umi_length + cprimer_start}" : "${cprimer_start}"

            PRESTO_MASKPRIMERS_EXTRACT_R1(
                ch_reads_R1,
                cprimer_start_var,
                primer_r1_extract_len,
                primer_r1_mask_mode,
                barcode_R1,
                "R1"
            )

            def barcode_R2 = (umi_position == "R2") ? true : false
            def vprimer_start_var = (umi_position == "R2") ? "${umi_length + vprimer_start}" : "${vprimer_start}"

            PRESTO_MASKPRIMERS_EXTRACT_R2(
                ch_reads_R2,
                vprimer_start_var,
                primer_r2_extract_len,
                primer_r2_mask_mode,
                barcode_R2,
                "R2"
            )
        } else if (cprimer_position == "R2") {
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false
            def vprimer_start_var = (index_file | umi_position == 'R1') ? "${umi_length + vprimer_start}" : "${vprimer_start}"

            PRESTO_MASKPRIMERS_EXTRACT_R1(
                ch_reads_R1,
                vprimer_start_var,
                primer_r1_extract_len,
                primer_r1_mask_mode,
                barcode_R1,
                "R1"
            )

            def barcode_R2 = (umi_position == "R2") ? true : false
            def cprimer_start_var = (umi_position == "R2") ? "${umi_length + cprimer_start}" : "${cprimer_start}"

            PRESTO_MASKPRIMERS_EXTRACT_R2(
                ch_reads_R2,
                cprimer_start_var,
                primer_r2_extract_len,
                primer_r2_mask_mode,
                barcode_R2,
                "R2"
            )
        } else {
            error "Error in determining cprimer position. Please choose R1 or R2."
        }


        // Merge again R1 and R2 by sample ID.
        ch_maskprimers_reads_R1 = PRESTO_MASKPRIMERS_EXTRACT_R1.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads_R2 = PRESTO_MASKPRIMERS_EXTRACT_R2.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads = ch_maskprimers_reads_R1.join(ch_maskprimers_reads_R2)
                                                        .map{ it -> [it[1], it[2], it[4]] }

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_EXTRACT_R1.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_EXTRACT_R2.out.logs)

        def barcode_position = (umi_position == "R1" | index_file) ? "R1" : "R2"
        PRESTO_PAIRSEQ_EXTRACT(
            ch_maskprimers_reads,
            barcode_position
        )
        ch_for_clustersets = PRESTO_PAIRSEQ_EXTRACT.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_EXTRACT.out.logs

    }  else {

        if (cprimer_position == "R1") {
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false
            def cprimer_start_var = (index_file | umi_position == 'R1') ? "${umi_length + cprimer_start}" : "${cprimer_start}"
            PRESTO_MASKPRIMERS_SCORE_UMI_R1(
                ch_reads_R1,
                ch_cprimers.collect(),
                cprimer_start_var,
                barcode_R1,
                primer_r1_maxerror,
                primer_r1_mask_mode,
                false,
                "R1"
            )

            def barcode_R2 = (umi_position == "R2") ? true : false
            def vprimer_start_var = (umi_position == "R2") ? "${umi_length + vprimer_start}" : "${vprimer_start}"
            PRESTO_MASKPRIMERS_SCORE_UMI_R2(
                ch_reads_R2,
                ch_vprimers.collect(),
                vprimer_start_var,
                barcode_R2,
                primer_r2_maxerror,
                primer_r2_mask_mode,
                primer_revpr,
                "R2"
            )
        } else if (cprimer_position == "R2") {
            def barcode_R1 = (index_file | umi_position == 'R1') ? true : false
            def vprimer_start_var = (index_file | umi_position == 'R1') ? "${umi_length + vprimer_start}" : "${vprimer_start}"
            PRESTO_MASKPRIMERS_SCORE_UMI_R1(
                ch_reads_R1,
                ch_vprimers.collect(),
                vprimer_start_var,
                barcode_R1,
                primer_r1_maxerror,
                primer_r1_mask_mode,
                primer_revpr,
                "R1"
            )
            def barcode_R2 = (umi_position == "R2") ? true : false
            def cprimer_start_var = (umi_position == "R2") ? "${umi_length + cprimer_start}" : "${cprimer_start}"

            PRESTO_MASKPRIMERS_SCORE_UMI_R2(
                ch_reads_R2,
                ch_cprimers.collect(),
                cprimer_start_var,
                barcode_R2,
                primer_r2_maxerror,
                primer_r2_mask_mode,
                primer_revpr,
                "R2"
            )
        } else {
            error "Error in determining cprimer position. Please choose R1 or R2."
        }


        // Merge again R1 and R2 by sample ID.
        ch_maskprimers_reads_R1 = PRESTO_MASKPRIMERS_SCORE_UMI_R1.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads_R2 = PRESTO_MASKPRIMERS_SCORE_UMI_R2.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}
        ch_maskprimers_reads = ch_maskprimers_reads_R1.join(ch_maskprimers_reads_R2)
                                                        .map{ it -> [it[1], it[2], it[4]] }

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_SCORE_UMI_R1.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_SCORE_UMI_R2.out.logs)

        // Pre-consensus pair
        def barcode_position = (umi_position == "R1" | index_file) ? "R1" : "R2"
        PRESTO_PAIRSEQ_UMI (
            ch_maskprimers_reads,
            barcode_position
        )

        ch_for_clustersets = PRESTO_PAIRSEQ_UMI.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_UMI.out.logs

    }

    if (cluster_sets) {

        // Cluster sequences by similarity
        PRESTO_CLUSTERSETS_UMI (
            ch_for_clustersets
        )

        // Annotate cluster into barcode field
        PRESTO_PARSE_CLUSTER_UMI (
            PRESTO_CLUSTERSETS_UMI.out.reads
        )
        ch_clustersets_logs = PRESTO_CLUSTERSETS_UMI.out.logs.collect()

        // Combining the split cluster+UMI combination
        PRESTO_PAIRSEQ_CLUSTERSETS (
            PRESTO_PARSE_CLUSTER_UMI.out.reads,
            "clustersets"
        )

        ch_for_buildconsensus = PRESTO_PAIRSEQ_CLUSTERSETS.out.reads

    } else {
        ch_for_buildconsensus = ch_for_clustersets
        ch_clustersets_logs = channel.empty()
    }

    // Build consensus of sequences with same UMI barcode
    if (maskprimers_align_race) {
        // Only consider C_PRIMER frequency when building consensus
        PRESTO_BUILDCONSENSUS_ALIGN_RACE (
            ch_for_buildconsensus,
            buildconsensus_maxerror,
            buildconsensus_maxgap,
            true, // use homogeneous C primer
            false, // do not use homogeneous linker
            primer_consensus,
            cluster_sets
        )
        ch_postconsensus = PRESTO_BUILDCONSENSUS_ALIGN_RACE.out.reads
        ch_readumi = PRESTO_BUILDCONSENSUS_ALIGN_RACE.out.log_tab
        ch_buildconsensus_logs = PRESTO_BUILDCONSENSUS_ALIGN_RACE.out.logs
    } else if (maskprimers_extract) {
        // Do not consider primers when building consensus
        PRESTO_BUILDCONSENSUS_EXTRACT(
            ch_for_buildconsensus,
            buildconsensus_maxerror,
            buildconsensus_maxgap,
            false,
            false,
            primer_consensus,
            cluster_sets
        )
        ch_postconsensus = PRESTO_BUILDCONSENSUS_EXTRACT.out.reads
        ch_readumi = PRESTO_BUILDCONSENSUS_EXTRACT.out.log_tab
        ch_buildconsensus_logs = PRESTO_BUILDCONSENSUS_EXTRACT.out.logs
    } else {
        // Consider both primers frequency when building consensus
        PRESTO_BUILDCONSENSUS_UMI (
            ch_for_buildconsensus,
            buildconsensus_maxerror,
            buildconsensus_maxgap,
            true,
            true,
            primer_consensus,
            cluster_sets
        )
        ch_postconsensus = PRESTO_BUILDCONSENSUS_UMI.out.reads
        ch_readumi = PRESTO_BUILDCONSENSUS_UMI.out.log_tab
        ch_buildconsensus_logs = PRESTO_BUILDCONSENSUS_UMI.out.logs
    }

    if (cluster_sets) {
        ch_readumi = PRESTO_CLUSTERSETS_UMI.out.log_tab
    }

    // Post-consensus pair
    PRESTO_POSTCONSENSUS_PAIRSEQ_UMI (
        ch_postconsensus
    )

    if (assemblepairs_sequential){
        // Assemble read pairs sequential
        PRESTO_ASSEMBLEPAIRS_SEQUENTIAL (
            PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.reads,
            ch_igblast.collect()
        )
        ch_assemblepairs_reads = PRESTO_ASSEMBLEPAIRS_SEQUENTIAL.out.reads
        ch_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_SEQUENTIAL.out.logs
    } else {
        // Assemble read pairs align
        PRESTO_ASSEMBLEPAIRS_UMI (
            PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.reads
        )
        ch_assemblepairs_reads = PRESTO_ASSEMBLEPAIRS_UMI.out.reads
        ch_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_UMI.out.logs
    }


    if (align_cregion) {
        PRESTO_ALIGN_CREGION(
            ch_assemblepairs_reads,
            ch_internal_cregion.collect(),
            cregion_maxlen,
            false,
            cregion_maxerror,
            cregion_mask_mode,
            false,
            "Cregion"
        )
        ch_parseheaders_reads = PRESTO_ALIGN_CREGION.out.reads
    } else {
        ch_parseheaders_reads = ch_assemblepairs_reads
    }

    // Generate QC stats after reads paired and filtered but before collapsed
    FASTQC_POSTASSEMBLY_UMI (
        ch_assemblepairs_reads
    )

    // Combine UMI duplicate count
    PRESTO_PARSEHEADERS_COLLAPSE_UMI (
        ch_parseheaders_reads
    )

    // Annotate primer fields and collapse duplicates
    if (maskprimers_align_race) {

        // Rename primer field to CREGION
        PRESTO_PARSEHEADERS_CREGION (
            PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads
        )

        // Collapse duplicates
        PRESTO_COLLAPSESEQ_ALIGN (
            PRESTO_PARSEHEADERS_CREGION.out.reads
        )
        ch_collapsed = PRESTO_COLLAPSESEQ_ALIGN.out.reads
        ch_collapse_logs = PRESTO_COLLAPSESEQ_ALIGN.out.logs

    } else if (maskprimers_extract) {

        // Do not annotate primers
        PRESTO_COLLAPSESEQ_UMI (
            PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads
        )
        ch_collapsed = PRESTO_COLLAPSESEQ_UMI.out.reads
        ch_collapse_logs = PRESTO_COLLAPSESEQ_UMI.out.logs

    } else {

        // Annotate primers in C_PRIMER and V_PRIMER field
        PRESTO_PARSEHEADERS_PRIMERS_UMI (
            PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads,
            cprimer_position
        )

        if (align_cregion) {
            PRESTO_COLLAPSESEQ_CREGION (
                PRESTO_PARSEHEADERS_PRIMERS_UMI.out.reads
            )
            ch_collapsed = PRESTO_COLLAPSESEQ_CREGION.out.reads
            ch_collapse_logs = PRESTO_COLLAPSESEQ_CREGION.out.logs
        } else {
            // Collapse duplicates
            PRESTO_COLLAPSESEQ_UMI (
                PRESTO_PARSEHEADERS_PRIMERS_UMI.out.reads
            )
            ch_collapsed = PRESTO_COLLAPSESEQ_UMI.out.reads
            ch_collapse_logs = PRESTO_COLLAPSESEQ_UMI.out.logs
        }
    }

    // Annotate metadata on read headers
    PRESTO_PARSEHEADERS_METADATA_UMI (
        ch_collapsed
    )

    // Filter out sequences with less than 2 representative duplicates with different UMIs
    PRESTO_SPLITSEQ_UMI (
        PRESTO_PARSEHEADERS_METADATA_UMI.out.reads
    )

    emit:
    fasta = PRESTO_SPLITSEQ_UMI.out.fasta
    versions = ch_versions
    fastp_reads_json = FASTP.out.json.collect{ meta,json -> json }
    fastp_reads_html = FASTP.out.html.collect{ meta,html -> html }
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_UMI.out.zip
    presto_filterseq_logs = PRESTO_FILTERSEQ_UMI.out.logs
    presto_maskprimers_logs = ch_maskprimers_logs.collect()
    presto_pairseq_logs = ch_pairseq_logs.collect()
    presto_clustersets_logs = ch_clustersets_logs
    presto_buildconsensus_logs = ch_buildconsensus_logs.collect()
    presto_postconsensus_pairseq_logs = PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.logs.collect()
    presto_assemblepairs_logs = ch_assemblepairs_logs.collect()
    presto_collapseseq_logs = ch_collapse_logs.collect()
    presto_splitseq_logs = PRESTO_SPLITSEQ_UMI.out.logs.collect()
    presto_UMIreads = ch_readumi.collect()
}
