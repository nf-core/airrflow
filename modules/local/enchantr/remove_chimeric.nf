process REMOVE_CHIMERIC {
    tag "$meta.id"
    label 'immcantation'
    label 'enchantr'
    label 'process_high'
    label 'process_long'

    // TODO: update container
    container "immcantation/suite:devel"


    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(imgt_base)

    output:
    tuple val(meta), path("*chimera-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" //, emit: chimera_report

    script:
    """
    Rscript -e "enchantr:::enchantr_report('chimera_analysis', report_params=list('input'='${tab}','outdir'=getwd(), 'nproc'=${task.cpus},'outname'='${meta.id}', 'log'='${meta.id}_chimeric_command_log'))"
    mv enchantr ${meta.id}_chimera_report
    """
}
