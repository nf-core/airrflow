// Include statements
include { MERGE_UMI                                      } from '../../modules/local/merge_UMI'
include { RENAME_FASTQ        as RENAME_FASTQ_UMI        } from '../../modules/local/rename_fastq'
include { GUNZIP              as GUNZIP_UMI              } from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_UMI } from '../../modules/local/fastqc_postassembly'
include { FASTP                                          } from '../../modules/nf-core/fastp/main'


//PRESTO
include { PRESTO_FILTERSEQ      as  PRESTO_FILTERSEQ_UMI }      from '../../modules/local/presto/presto_filterseq'
include { PRESTO_MASKPRIMERS    as  PRESTO_MASKPRIMERS_UMI }    from '../../modules/local/presto/presto_maskprimers'
include { PRESTO_MASKPRIMERS_ALIGN }                            from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_MASKPRIMERS_EXTRACT }                          from '../../modules/local/presto/presto_maskprimers_extract'
include { PRESTO_MASKPRIMERS_ALIGN as PRESTO_ALIGN_CREGION }    from '../../modules/local/presto/presto_maskprimers_align'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_UMI }        from '../../modules/local/presto/presto_pairseq'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_ALIGN } from '../../modules/local/presto/presto_pairseq'
include { PRESTO_CLUSTERSETS    as  PRESTO_CLUSTERSETS_UMI }    from '../../modules/local/presto/presto_clustersets'
include { PRESTO_PARSE_CLUSTER  as  PRESTO_PARSE_CLUSTER_UMI }  from '../../modules/local/presto/presto_parse_cluster'
include { PRESTO_BUILDCONSENSUS as  PRESTO_BUILDCONSENSUS_UMI } from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_BUILDCONSENSUS as PRESTO_BUILDCONSENSUS_ALIGN } from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_POSTCONSENSUS_PAIRSEQ as PRESTO_POSTCONSENSUS_PAIRSEQ_UMI }    from '../../modules/local/presto/presto_postconsensus_pairseq'
include { PRESTO_ASSEMBLEPAIRS  as  PRESTO_ASSEMBLEPAIRS_UMI }  from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_ASSEMBLEPAIRS_SEQUENTIAL }                     from '../../modules/local/presto/presto_assemblepairs_sequential'
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_COLLAPSE_UMI } from '../../modules/local/presto/presto_parseheaders'
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_CREGION }   from '../../modules/local/presto/presto_parseheaders'
include { PRESTO_PARSEHEADERS_PRIMERS   as PRESTO_PARSEHEADERS_PRIMERS_UMI }    from '../../modules/local/presto/presto_parseheaders_primers'
include { PRESTO_PARSEHEADERS_METADATA  as PRESTO_PARSEHEADERS_METADATA_UMI }   from '../../modules/local/presto/presto_parseheaders_metadata'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_UMI }     from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_ALIGN }   from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_CREGION } from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_SPLITSEQ       as PRESTO_SPLITSEQ_UMI}         from '../../modules/local/presto/presto_splitseq'


workflow PRESTO_UMI {
    take:
    ch_reads       // channel: [ val(meta), [ reads ] ]
    ch_cprimers    // channel: [ cprimers.fasta ]
    ch_vprimers    // channel: [ vprimers.fasta ]
    ch_adapter_fasta // channel: [ adapters.fasta ]
    ch_internal_cregion // channel: [ internal_cregions.fasta ]
    ch_igblast

    main:

    ch_versions = Channel.empty()

    // Merge UMI from index file to R1 if provided
    if (params.index_file) {

        // ch for fastp reads R1 R2
        ch_reads.map{ meta, reads -> [meta, [reads[0], reads[1]]] }
                .set{ ch_reads_R1_R2 }

        // Fastp reads R1 R2
        save_merged = false
        FASTP (
            ch_reads_R1_R2,
            ch_adapter_fasta,
            params.save_trimmed,
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

        MERGE_UMI ( ch_meta_R1_R2_index )
        ch_gunzip = MERGE_UMI.out.reads
        ch_versions = ch_versions.mix(MERGE_UMI.out.versions)


    } else {

        // Fastp reads
        save_merged = false
        FASTP (
            ch_reads,
            ch_adapter_fasta,
            params.save_trimmed,
            save_merged
        )
        ch_versions = ch_versions.mix(FASTP.out.versions)

        ch_rename_fastq_umi = FASTP.out.reads.map{ meta,reads -> [meta, reads[0], reads[1]] }

        RENAME_FASTQ_UMI ( ch_rename_fastq_umi )
        ch_gunzip = RENAME_FASTQ_UMI.out.reads

    }

    // gunzip fastq.gz to fastq
    GUNZIP_UMI ( ch_gunzip )
    ch_versions = ch_versions.mix(GUNZIP_UMI.out.versions)

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_UMI ( GUNZIP_UMI.out.reads )
    ch_versions = ch_versions.mix(PRESTO_FILTERSEQ_UMI.out.versions)

    // Mask primers
    if (params.maskprimers_align) {

        ch_reads_R1 = PRESTO_FILTERSEQ_UMI.out.reads
                                            .map{ reads -> [reads[0], reads[1]] }.dump(tag: 'ch_reads_R1')
        ch_reads_R2 = PRESTO_FILTERSEQ_UMI.out.reads
                                            .map{ reads -> [reads[0], reads[2]] }.dump(tag: 'ch_reads_R2')
        PRESTO_MASKPRIMERS_ALIGN(
            ch_reads_R1,
            ch_cprimers.collect(),
            params.primer_maxlen,
            params.primer_r1_maxerror,
            params.primer_mask_mode
        )
        PRESTO_MASKPRIMERS_EXTRACT(
            ch_reads_R2
        )

        ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_ALIGN.out.versions)
        // Merge again R1 and R2 by sample ID.
        ch_maskprimers_reads_R1 = PRESTO_MASKPRIMERS_ALIGN.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}.dump(tag: 'ch_maskprimers_reads_R1')
        ch_maskprimers_reads_R2 = PRESTO_MASKPRIMERS_EXTRACT.out.reads.map{ reads -> [reads[0].id, reads[0], reads[1]]}.dump(tag: 'ch_maskprimers_reads_R2')
        ch_maskprimers_reads = ch_maskprimers_reads_R1.join(ch_maskprimers_reads_R2)
                                                        .map{ it -> [it[1], it[2], it[4]] }.dump(tag: 'ch_maskprimers_reads_after_remerge')

        ch_maskprimers_logs = PRESTO_MASKPRIMERS_ALIGN.out.logs
        ch_maskprimers_logs = ch_maskprimers_logs.mix(PRESTO_MASKPRIMERS_EXTRACT.out.logs)

        PRESTO_PAIRSEQ_ALIGN( ch_maskprimers_reads )
        ch_versions  = ch_versions.mix(PRESTO_PAIRSEQ_ALIGN.out.versions)
        ch_for_clustersets = PRESTO_PAIRSEQ_ALIGN.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_ALIGN.out.logs

    } else {

        PRESTO_MASKPRIMERS_UMI (
            PRESTO_FILTERSEQ_UMI.out.reads,
            ch_cprimers.collect(),
            ch_vprimers.collect()
        )
        ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_UMI.out.versions)
        ch_maskprimers_logs = PRESTO_MASKPRIMERS_UMI.out.logs

        // Pre-consensus pair
        PRESTO_PAIRSEQ_UMI (
            PRESTO_MASKPRIMERS_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_PAIRSEQ_UMI.out.versions)
        ch_for_clustersets = PRESTO_PAIRSEQ_UMI.out.reads
        ch_pairseq_logs = PRESTO_PAIRSEQ_UMI.out.logs

    }

    if (params.cluster_sets) {

        // Cluster sequences by similarity
        PRESTO_CLUSTERSETS_UMI (
            ch_for_clustersets
        )
        ch_versions = ch_versions.mix(PRESTO_CLUSTERSETS_UMI.out.versions)

        // Annotate cluster into barcode field
        PRESTO_PARSE_CLUSTER_UMI (
            PRESTO_CLUSTERSETS_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_PARSE_CLUSTER_UMI.out.versions)
        ch_for_buildconsensus = PRESTO_PARSE_CLUSTER_UMI.out.reads
        ch_clustersets_logs = PRESTO_CLUSTERSETS_UMI.out.logs.collect()

    } else {
        ch_for_buildconsensus = ch_for_clustersets
        ch_clustersets_logs = Channel.empty()
    }

    // Build consensus of sequences with same UMI barcode
    if (params.maskprimers_align) {
        PRESTO_BUILDCONSENSUS_ALIGN (
            ch_for_buildconsensus
        )
        ch_versions = ch_versions.mix(PRESTO_BUILDCONSENSUS_ALIGN.out.versions)
        ch_postconsensus = PRESTO_BUILDCONSENSUS_ALIGN.out.reads
        ch_buildconsensus_logs = PRESTO_BUILDCONSENSUS_ALIGN.out.logs
    } else {
        PRESTO_BUILDCONSENSUS_UMI (
            ch_for_buildconsensus
        )
        ch_versions = ch_versions.mix(PRESTO_BUILDCONSENSUS_UMI.out.versions)
        ch_postconsensus = PRESTO_BUILDCONSENSUS_UMI.out.reads
        ch_buildconsensus_logs = PRESTO_BUILDCONSENSUS_UMI.out.logs
    }

    // Post-consensus pair
    PRESTO_POSTCONSENSUS_PAIRSEQ_UMI (
        ch_postconsensus
    )
    ch_versions = ch_versions.mix(PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.versions)

    if (params.assemblepairs_sequential){
        // Assemble read pairs sequential
        PRESTO_ASSEMBLEPAIRS_SEQUENTIAL (
            PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.reads,
            ch_igblast.collect()
        )
        ch_versions = ch_versions.mix(PRESTO_ASSEMBLEPAIRS_SEQUENTIAL.out.versions)
        ch_assemblepairs_reads = PRESTO_ASSEMBLEPAIRS_SEQUENTIAL.out.reads
        ch_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_SEQUENTIAL.out.logs
    } else {
        // Assemble read pairs align
        PRESTO_ASSEMBLEPAIRS_UMI (
            PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_ASSEMBLEPAIRS_UMI.out.versions)
        ch_assemblepairs_reads = PRESTO_ASSEMBLEPAIRS_UMI.out.reads
        ch_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_UMI.out.logs
    }


    if (params.align_cregion) {
        PRESTO_ALIGN_CREGION(
            ch_assemblepairs_reads,
            ch_internal_cregion.collect(),
            params.cregion_maxlen,
            params.cregion_maxerror,
            params.cregion_mask_mode
        )
        ch_parseheaders_reads = PRESTO_ALIGN_CREGION.out.reads
    } else {
        ch_parseheaders_reads = ch_assemblepairs_reads
    }

    // Generate QC stats after reads paired and filtered but before collapsed
    FASTQC_POSTASSEMBLY_UMI (
        ch_assemblepairs_reads
    )
    ch_versions = ch_versions.mix(FASTQC_POSTASSEMBLY_UMI.out.versions)

    // Combine UMI duplicate count
    PRESTO_PARSEHEADERS_COLLAPSE_UMI (
        ch_parseheaders_reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.versions)

    // Annotate primer fields and collapse duplicates
    if (params.maskprimers_align) {
        // Rename primer field to CREGION
        PRESTO_PARSEHEADERS_CREGION (
            PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_CREGION.out.versions)

        // Collapse duplicates
        PRESTO_COLLAPSESEQ_ALIGN (
            PRESTO_PARSEHEADERS_CREGION.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_COLLAPSESEQ_ALIGN.out.versions)
        ch_collapsed = PRESTO_COLLAPSESEQ_ALIGN.out.reads
        ch_collapse_logs = PRESTO_COLLAPSESEQ_ALIGN.out.logs

    } else {
        // Annotate primers in C_PRIMER and V_PRIMER field
        PRESTO_PARSEHEADERS_PRIMERS_UMI (
            PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_PRIMERS_UMI.out.versions)

        if (params.align_cregion) {
            PRESTO_COLLAPSESEQ_CREGION (
                PRESTO_PARSEHEADERS_PRIMERS_UMI.out.reads
            )
            ch_versions = ch_versions.mix(PRESTO_COLLAPSESEQ_CREGION.out.versions)
            ch_collapsed = PRESTO_COLLAPSESEQ_CREGION.out.reads
            ch_collapse_logs = PRESTO_COLLAPSESEQ_CREGION.out.logs
        } else {
            // Collapse duplicates
            PRESTO_COLLAPSESEQ_UMI (
                PRESTO_PARSEHEADERS_PRIMERS_UMI.out.reads
            )
            ch_versions = ch_versions.mix(PRESTO_COLLAPSESEQ_UMI.out.versions)
            ch_collapsed = PRESTO_COLLAPSESEQ_UMI.out.reads
            ch_collapse_logs = PRESTO_COLLAPSESEQ_UMI.out.logs
        }
    }

    // Annotate metadata on read headers
    PRESTO_PARSEHEADERS_METADATA_UMI (
        ch_collapsed
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_METADATA_UMI.out.versions)

    // Filter out sequences with less than 2 representative duplicates with different UMIs
    PRESTO_SPLITSEQ_UMI (
        PRESTO_PARSEHEADERS_METADATA_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_SPLITSEQ_UMI.out.versions)

    emit:
    fasta = PRESTO_SPLITSEQ_UMI.out.fasta
    software = ch_versions
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
}
