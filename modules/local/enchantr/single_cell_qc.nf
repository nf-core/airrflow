process SINGLE_CELL_QC {
    tag 'all_single_cell'
    label 'immcantation'
    label 'enchantr'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    path tabs

    output:
    path("*scqc-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path("*_report"), emit: report
    path("versions.yml"), emit: versions

    script:
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('single_cell_qc', report_params=list('input'='tabs.txt','outdir'=getwd(), 'outname'='all_reps', 'log'='all_reps_scqc_command_log'))"
    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    mv enchantr all_reps_scqc_report
    """
}
