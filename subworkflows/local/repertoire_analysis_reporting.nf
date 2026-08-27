include { PARSE_LOGS } from '../../modules/local/parse_logs'
include { REPORT_FILE_SIZE } from '../../modules/local/enchantr/report_file_size'
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
    ch_input_check_logs
    ch_reassign_logs
    ch_changeo_makedb_logs
    ch_vdj_annotation_logs
    ch_bulk_qc_and_filter_logs
    ch_sc_qc_and_filter_logs
    ch_input // Input samplesheet
    ch_report_rmd // Report Rmarkdown file
    ch_report_css // Report CSS file
    ch_report_logo // Logo to be displayed in report
    ch_metadata // Validated samplesheet
    mode
    library_generation_method
    umi_length
    cluster_sets

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

    REPORT_FILE_SIZE(
        ch_logs.collect().ifEmpty([]),
        ch_metadata,
        ch_logs_tabs
    )

    AIRRFLOW_REPORT(
        ch_parsed_logs.collect().ifEmpty([]),
        REPORT_FILE_SIZE.out.table.collect().ifEmpty([]),
        ch_report_rmd,
        ch_report_css,
        ch_report_logo
    )
}
