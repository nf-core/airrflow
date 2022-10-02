process COLLAPSE_DUPLICATES {
    tag "all_bulk_reps"
    label 'immcantation'
    label 'enchantr'
    label 'process_long'

    cache  'lenient'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    path(tabs) // tuple [val(meta), sequence tsv in AIRR format ]

    output:
    tuple val(meta), path("*collapse-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" // , emit: duplicates_report
    path "versions.yml" , emit: versions

    script:
    meta=[]
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('collapse_duplicates', report_params=list('input'='tabs.txt','collapseby'='${params.collapseby}','outdir'=getwd(), 'nproc'=${task.cpus},'outname'='all_reps', 'log'='all_reps_collapse_command_log'))"
    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    mv enchantr all_reps_collapse_report
    """
}
