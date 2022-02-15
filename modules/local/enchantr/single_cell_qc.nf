process SINGLE_CELL_QC {
    tag "repertoire_all"
    label 'immcantation'
    label 'enchantr'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    //tuple val(meta), path(tab) // sequence tsv in AIRR format
    path tabs

    output:
    tuple val(meta), path("*scqc-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" //, em`it: chimera_report

    script:
    meta=[]
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('single_cell_qc', report_params=list('input'='tabs.txt','outdir'=getwd(), 'outname'='all_reps', 'log'='all_reps_scqc_command_log'))"
    mv enchantr all_reps_scqc_report
    """
}
