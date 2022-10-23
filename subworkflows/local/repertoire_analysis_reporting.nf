include { PARSE_LOGS } from '../../modules/local/parse_logs.nf'
include { ALAKAZAM_SHAZAM_REPERTOIRES  } from '../../modules/local/alakazam/alakazam_shazam_repertoires'

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
    ch_changeo_makedb_logs
    ch_repertoires
    ch_input
    ch_report_rmd
    ch_report_css
    ch_report_logo

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


    ALAKAZAM_SHAZAM_REPERTOIRES(
        ch_repertoires,
        ch_parsed_logs.collect(),
        ch_report_rmd,
        ch_report_css,
        ch_report_logo
    )
    ch_versions = ch_versions.mix(PARSE_LOGS.out.versions)

    emit:
    versions = ch_versions
}
