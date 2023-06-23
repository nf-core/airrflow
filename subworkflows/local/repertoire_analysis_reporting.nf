include { PARSE_LOGS } from '../../modules/local/parse_logs.nf'
include { REPORT_FILE_SIZE } from '../../modules/local/enchantr/report_file_size.nf'
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
    ch_reassign_logs
    ch_changeo_makedb_logs
    ch_vdj_annotation_logs
    ch_bulk_qc_and_filter_logs
    ch_sc_qc_and_filter_logs
    ch_clonal_analysis_logs
    ch_repertoires // Repertoire tsv files from clonal analysis process
    ch_input // Input samplesheet
    ch_report_rmd // Report Rmarkdown file
    ch_report_css // Report CSS file
    ch_report_logo // Logo to be displayed in report
    ch_metadata // Validated samplesheet

    main:
    ch_versions = Channel.empty()

    if (params.mode == "fastq") {
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
            ch_input
        )
        ch_versions = ch_versions.mix(PARSE_LOGS.out.versions)
        ch_parsed_logs = PARSE_LOGS.out.logs

    } else {
        ch_parsed_logs = Channel.empty()
    }

    ch_logs = ch_vdj_annotation_logs.mix(ch_bulk_qc_and_filter_logs,
                                        ch_reassign_logs,
                                        ch_sc_qc_and_filter_logs,
                                        ch_clonal_analysis_logs)
    ch_logs_tabs =  ch_logs.collect()
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .dump(tag: 'ch_logs_tabs')
                        .collectFile(name: 'all_logs_tabs.txt', newLine: true)

    REPORT_FILE_SIZE(
        ch_logs.collect().ifEmpty([]),
        ch_metadata,
        ch_logs_tabs
    )
    ch_versions = ch_versions.mix(REPORT_FILE_SIZE.out.versions)

    AIRRFLOW_REPORT(
        ch_repertoires,
        ch_parsed_logs.collect().ifEmpty([]),
        REPORT_FILE_SIZE.out.table.ifEmpty([]),
        ch_report_rmd,
        ch_report_css,
        ch_report_logo
    )
    ch_versions = ch_versions.mix(AIRRFLOW_REPORT.out.versions)

    emit:
    versions = ch_versions
}
