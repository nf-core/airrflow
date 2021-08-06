def modules = params.modules.clone()

include { GUNZIP            as GUNZIP_SANS_UMI }         from '../../modules/local/gunzip'                    addParams( options: [:] )
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_SANS_UMI } from '../../modules/local/fastqc_postassembly'                                addParams( options: [:] )

//PRESTO
include { PRESTO_ASSEMBLEPAIRS  as  PRESTO_ASSEMBLEPAIRS_SANS_UMI }  from '../../modules/local/presto/presto_assemblepairs'                         addParams( options: modules['presto_assemblepairs_sans_umi'] )
include { PRESTO_FILTERSEQ_POSTASSEMBLY      as  PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI }      from '../../modules/local/presto/presto_filterseq_postassembly'                             addParams( options: modules['presto_filterseq_postassembly_sans_umi'] )
include { PRESTO_MASKPRIMERS_POSTASSEMBLY    as  PRESTO_MASKPRIMERS_POSTASSEMBLY_SANS_UMI }    from '../../modules/local/presto/presto_maskprimers_postassembly'                           addParams( options: modules['presto_maskprimers_postassembly_sans_umi'] )
include { PRESTO_PARSEHEADERS   as  PRESTO_PARSEHEADERS_COLLAPSE_SANS_UMI } from '../../modules/local/presto/presto_parseheaders'                   addParams( options: modules['presto_parseheaders_collapse_sans_umi'] )
include { PRESTO_PARSEHEADERS_PRIMERS   as PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI }    from '../../modules/local/presto/presto_parseheaders_primers'      addParams( options: modules['presto_parseheaders_primers_sans_umi'] )
include { PRESTO_PARSEHEADERS_METADATA  as PRESTO_PARSEHEADERS_METADATA_SANS_UMI }   from '../../modules/local/presto/presto_parseheaders_metadata'     addParams( options: modules['presto_parseheaders_metadata'] )
include { PRESTO_COLLAPSESEQ    as PRESTO_COLLAPSESEQ_SANS_UMI }     from '../../modules/local/presto/presto_collapseseq'                           addParams( options: modules['presto_collapseseq_sans_umi'] )
include { PRESTO_SPLITSEQ       as PRESTO_SPLITSEQ_SANS_UMI}         from '../../modules/local/presto/presto_splitseq'                              addParams( options: modules['presto_splitseq_sans_umi'] )


workflow PRESTO_SANS_UMI {
    take:
    ch_reads       // channel: [ val(meta), [ reads ] ]
    ch_cprimers    // channel: [ cprimers.fasta ]
    ch_vprimers    // channel: [ vprimers.fasta ]

    main:

    ch_software_versions = Channel.empty()
    ch_gunzip = ch_reads

    // gunzip fastq.gz to fastq
    GUNZIP_SANS_UMI ( ch_gunzip )
    ch_software_versions = ch_software_versions.mix(GUNZIP_SANS_UMI.out.version.first().ifEmpty(null))

    // Assemble read pairs
    PRESTO_ASSEMBLEPAIRS_SANS_UMI (
        GUNZIP_SANS_UMI.out.reads
    )

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI (
        PRESTO_ASSEMBLEPAIRS_SANS_UMI.out.reads
    )
    ch_software_versions = ch_software_versions.mix(PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.version.first().ifEmpty(null))

    // Mask primers
    PRESTO_MASKPRIMERS_POSTASSEMBLY_SANS_UMI (
        PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.reads,
        ch_cprimers.collect(),
        ch_vprimers.collect()
    )

    // Generate QC stats after reads paired and filtered but before collapsed
    FASTQC_POSTASSEMBLY_SANS_UMI (
        PRESTO_MASKPRIMERS_POSTASSEMBLY_SANS_UMI.out.reads
    )

    // Collapse duplicates into counts
    PRESTO_PARSEHEADERS_COLLAPSE_SANS_UMI (
        PRESTO_MASKPRIMERS_POSTASSEMBLY_SANS_UMI.out.reads
    )

    // Annotate primers in C_PRIMER and V_PRIMER field
    PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI (
        PRESTO_PARSEHEADERS_COLLAPSE_SANS_UMI.out.reads
    )

    // Annotate metadata on primer headers
    PRESTO_PARSEHEADERS_METADATA_SANS_UMI (
        PRESTO_PARSEHEADERS_PRIMERS_SANS_UMI.out.reads
    )

    // Mark and count duplicate sequences with different UMI barcodes (DUPCOUNT)
    PRESTO_COLLAPSESEQ_SANS_UMI (
        PRESTO_PARSEHEADERS_METADATA_SANS_UMI.out.reads
    )

    // Filter out sequences with less than 2 representative duplicates with different UMIs
    PRESTO_SPLITSEQ_SANS_UMI (
        PRESTO_COLLAPSESEQ_SANS_UMI.out.reads
    )

    emit:
    fasta = PRESTO_SPLITSEQ_SANS_UMI.out.fasta
    software = ch_software_versions
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_SANS_UMI.out.zip
    presto_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_SANS_UMI.out.logs.collect()
    presto_filterseq_logs = PRESTO_FILTERSEQ_POSTASSEMBLY_SANS_UMI.out.logs
    presto_maskprimers_logs = PRESTO_MASKPRIMERS_POSTASSEMBLY_SANS_UMI.out.logs.collect()
    presto_collapseseq_logs = PRESTO_COLLAPSESEQ_SANS_UMI.out.logs.collect()
    presto_splitseq_logs = PRESTO_SPLITSEQ_SANS_UMI.out.logs.collect()
}
