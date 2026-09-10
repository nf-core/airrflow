include { PARSE_LOGS } from '../../modules/local/parse_logs'
include { REPORT_FILE_SIZE } from '../../modules/local/enchantr/report_file_size'
include { GET_READS_PER_UMI } from '../../modules/local/get_reads_per_umi'
include { AIRRFLOW_REPORT  } from '../../modules/local/airrflow_report/airrflow_report'

workflow REPERTOIRE_ANALYSIS_REPORTING {

    take:
    ch_presto_filterseq_logs
    ch_presto_maskprimers_logs
    ch_presto_pairseq_logs
    ch_presto_clustersets_logs
    ch_presto_buildconsensus_logs
    ch_presto_postconsensus_pairseq_logs
    ch_presto_assemblepairs_logs
    ch_presto_collapseseq_logs
    ch_presto_splitseq_logs
    ch_presto_UMIreads
    ch_fasta
    ch_input_check_logs
    ch_reassign_logs
    ch_changeo_makedb_logs
    ch_vdj_annotation_logs
    ch_bulk_qc_and_filter_logs
    ch_sc_qc_and_filter_logs
    ch_repertoires // Repertoire tsv files from clonal analysis process
    ch_input // Input samplesheet
    ch_report_rmd // Report Rmarkdown file
    ch_report_css // Report CSS file
    ch_report_logo // Logo to be displayed in report
    ch_metadata // Validated samplesheet
    mode
    library_generation_method
    umi_length
    cluster_sets
    trust4_umi_read

    main:

    if (mode == "fastq" && library_generation_method != "sc_10x_genomics" && library_generation_method != "trust4" ) {
        PARSE_LOGS(
            ch_presto_filterseq_logs,
            ch_presto_maskprimers_logs,
            ch_presto_pairseq_logs,
            ch_presto_clustersets_logs,
            ch_presto_buildconsensus_logs,
            ch_presto_postconsensus_pairseq_logs,
            ch_presto_assemblepairs_logs,
            ch_presto_collapseseq_logs,
            ch_presto_splitseq_logs,
            ch_changeo_makedb_logs,
            ch_input,
            umi_length,
            cluster_sets
        )
        ch_parsed_logs = PARSE_LOGS.out.logs

    } else {
        ch_parsed_logs = channel.empty()
    }

    ch_logs = ch_vdj_annotation_logs.mix(
                                        ch_input_check_logs,
                                        ch_bulk_qc_and_filter_logs,
                                        ch_reassign_logs,
                                        ch_sc_qc_and_filter_logs)
    ch_logs_tabs =  ch_logs.collect()
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .collectFile(name: 'all_logs_tabs.txt', newLine: true)

    // if UMI protocol
    if ((library_generation_method == 'specific_pcr_umi') || (library_generation_method == 'dt_5p_race_umi')) {
        GET_READS_PER_UMI(
            ch_presto_UMIreads,
            [],
            library_generation_method
        )
        ch_reads_umi = GET_READS_PER_UMI.out.readUMI.collect()
    } else if ((library_generation_method == 'trust4')&&(trust4_umi_read ? trust4_umi_read : [])) {
        GET_READS_PER_UMI(
            [],
            ch_fasta,
            library_generation_method
        )
        ch_reads_umi = GET_READS_PER_UMI.out.readUMI.collect()
    } else {
        ch_reads_umi = []
    }

    REPORT_FILE_SIZE(
        ch_logs.collect().ifEmpty([]),
        ch_metadata,
        ch_logs_tabs
    )

    AIRRFLOW_REPORT(
        ch_repertoires,
        ch_parsed_logs.collect().ifEmpty([]),
        REPORT_FILE_SIZE.out.table.collect().ifEmpty([]),
        ch_report_rmd,
        ch_report_css,
        ch_report_logo,
        ch_reads_umi
    )
}
