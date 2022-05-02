// Include statements
include { MERGE_UMI }                               from '../../modules/local/merge_UMI'
include { RENAME_FASTQ      as RENAME_FASTQ_UMI }   from '../../modules/local/rename_fastq'
include { GUNZIP            as GUNZIP_UMI }         from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_UMI } from '../../modules/local/fastqc_postassembly'

//PRESTO
include { PRESTO_FILTERSEQ      as  PRESTO_FILTERSEQ_UMI }      from '../../modules/local/presto/presto_filterseq'
include { PRESTO_MASKPRIMERS    as  PRESTO_MASKPRIMERS_UMI }    from '../../modules/local/presto/presto_maskprimers'
include { PRESTO_PAIRSEQ        as  PRESTO_PAIRSEQ_UMI }        from '../../modules/local/presto/presto_pairseq'
include { PRESTO_CLUSTERSETS    as  PRESTO_CLUSTERSETS_UMI }    from '../../modules/local/presto/presto_clustersets'
include { PRESTO_PARSE_CLUSTER  as  PRESTO_PARSE_CLUSTER_UMI }  from '../../modules/local/presto/presto_parse_cluster'
include { PRESTO_BUILDCONSENSUS as  PRESTO_BUILDCONSENSUS_UMI } from '../../modules/local/presto/presto_buildconsensus'
include { PRESTO_POSTCONSENSUS_PAIRSEQ as PRESTO_POSTCONSENSUS_PAIRSEQ_UMI }    from '../../modules/local/presto/presto_postconsensus_pairseq'
include { PRESTO_ASSEMBLEPAIRS  as  PRESTO_ASSEMBLEPAIRS_UMI }  from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_COLLAPSE_UMI } from '../../modules/local/presto/presto_parseheaders'
include { PRESTO_PARSEHEADERS_PRIMERS   as PRESTO_PARSEHEADERS_PRIMERS_UMI }    from '../../modules/local/presto/presto_parseheaders_primers'
include { PRESTO_PARSEHEADERS_METADATA  as PRESTO_PARSEHEADERS_METADATA_UMI }   from '../../modules/local/presto/presto_parseheaders_metadata'
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_UMI }     from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_SPLITSEQ       as PRESTO_SPLITSEQ_UMI}         from '../../modules/local/presto/presto_splitseq'


workflow PRESTO_UMI {
    take:
    ch_reads       // channel: [ val(meta), [ reads ] ]
    ch_cprimers    // channel: [ cprimers.fasta ]
    ch_vprimers    // channel: [ vprimers.fasta ]

    main:

    ch_versions = Channel.empty()
    // Merge UMI from index file to R1 if provided
    if (params.index_file) {
        MERGE_UMI ( ch_reads )
        ch_gunzip = MERGE_UMI.out.reads
        ch_versions = ch_versions.mix(MERGE_UMI.out.versions.ifEmpty(null))
    } else {
        RENAME_FASTQ_UMI ( ch_reads )
        ch_gunzip = RENAME_FASTQ_UMI.out.reads
    }

    // gunzip fastq.gz to fastq
    GUNZIP_UMI ( ch_gunzip )
    ch_versions = ch_versions.mix(GUNZIP_UMI.out.versions.ifEmpty(null))

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_UMI ( GUNZIP_UMI.out.reads )
    ch_versions = ch_versions.mix(PRESTO_FILTERSEQ_UMI.out.versions.ifEmpty(null))

    // Mask primers
    PRESTO_MASKPRIMERS_UMI (
        PRESTO_FILTERSEQ_UMI.out.reads,
        ch_cprimers.collect(),
        ch_vprimers.collect()
    )
    ch_versions = ch_versions.mix(PRESTO_MASKPRIMERS_UMI.out.versions.ifEmpty(null))

    // Pre-consensus pair
    PRESTO_PAIRSEQ_UMI (
        PRESTO_MASKPRIMERS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_PAIRSEQ_UMI.out.versions.ifEmpty(null))

    if (params.cluster_sets) {

        // Cluster sequences by similarity
        PRESTO_CLUSTERSETS_UMI (
            PRESTO_PAIRSEQ_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_CLUSTERSETS_UMI.out.versions.ifEmpty(null))

        // Annotate cluster into barcode field
        PRESTO_PARSE_CLUSTER_UMI (
            PRESTO_CLUSTERSETS_UMI.out.reads
        )
        ch_versions = ch_versions.mix(PRESTO_PARSE_CLUSTER_UMI.out.versions.ifEmpty(null))
        ch_for_buildconsensus = PRESTO_PARSE_CLUSTER_UMI.out.reads
        ch_clustersets_logs = PRESTO_CLUSTERSETS_UMI.out.logs.collect()

    } else {
        ch_for_buildconsensus = PRESTO_PAIRSEQ_UMI.out.reads
        ch_clustersets_logs = Channel.empty()
    }

    // Build consensus of sequences with same UMI barcode
    PRESTO_BUILDCONSENSUS_UMI (
        ch_for_buildconsensus
    )
    ch_versions = ch_versions.mix(PRESTO_BUILDCONSENSUS_UMI.out.versions.ifEmpty(null))

    // Post-consensus pair
    PRESTO_POSTCONSENSUS_PAIRSEQ_UMI (
        PRESTO_BUILDCONSENSUS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.versions.ifEmpty(null))

    // Assemble read pairs
    PRESTO_ASSEMBLEPAIRS_UMI (
        PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_ASSEMBLEPAIRS_UMI.out.versions.ifEmpty(null))

    // Generate QC stats after reads paired and filtered but before collapsed
    FASTQC_POSTASSEMBLY_UMI (
        PRESTO_ASSEMBLEPAIRS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_POSTASSEMBLY_UMI.out.versions.ifEmpty(null))

    // Combine UMI duplicate count
    PRESTO_PARSEHEADERS_COLLAPSE_UMI (
        PRESTO_ASSEMBLEPAIRS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.versions.ifEmpty(null))

    // Annotate primers in C_PRIMER and V_PRIMER field
    PRESTO_PARSEHEADERS_PRIMERS_UMI (
        PRESTO_PARSEHEADERS_COLLAPSE_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_PRIMERS_UMI.out.versions.ifEmpty(null))

    // Annotate metadata on primer headers
    PRESTO_PARSEHEADERS_METADATA_UMI (
        PRESTO_PARSEHEADERS_PRIMERS_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_PARSEHEADERS_METADATA_UMI.out.versions.ifEmpty(null))

    // Mark and count duplicate sequences with different UMI barcodes (DUPCOUNT)
    PRESTO_COLLAPSESEQ_UMI (
        PRESTO_PARSEHEADERS_METADATA_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_COLLAPSESEQ_UMI.out.versions.ifEmpty(null))

    // Filter out sequences with less than 2 representative duplicates with different UMIs
    PRESTO_SPLITSEQ_UMI (
        PRESTO_COLLAPSESEQ_UMI.out.reads
    )
    ch_versions = ch_versions.mix(PRESTO_SPLITSEQ_UMI.out.versions.ifEmpty(null))

    emit:
    fasta = PRESTO_SPLITSEQ_UMI.out.fasta
    software = ch_versions
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_UMI.out.zip
    presto_filterseq_logs = PRESTO_FILTERSEQ_UMI.out.logs
    presto_maskprimers_logs = PRESTO_MASKPRIMERS_UMI.out.logs.collect()
    presto_pairseq_logs = PRESTO_PAIRSEQ_UMI.out.logs.collect()
    presto_clustersets_logs = ch_clustersets_logs
    presto_buildconsensus_logs = PRESTO_BUILDCONSENSUS_UMI.out.logs.collect()
    presto_postconsensus_pairseq_logs = PRESTO_POSTCONSENSUS_PAIRSEQ_UMI.out.logs.collect()
    presto_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_UMI.out.logs.collect()
    presto_collapseseq_logs = PRESTO_COLLAPSESEQ_UMI.out.logs.collect()
    presto_splitseq_logs = PRESTO_SPLITSEQ_UMI.out.logs.collect()
}
