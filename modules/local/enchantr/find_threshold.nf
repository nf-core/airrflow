process FIND_THRESHOLD {
    tag "all_reps"
    label 'immcantation'
    label 'enchantr'
    label 'process_high'
    label 'process_long'

    cache 'lenient'

    // TODO: update container
    container "immcantation/suite:devel"


    input:
    path tab // sequence tsv in AIRR format

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "*_threshold-summary.tsv"
    path "*_threshold-mean.tsv", emit: mean_threshold

    script:
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', report_params=list('input'='${tab.join(',')}','cloneby'='${params.cloneby}','singlecell'='${params.singlecell}','outdir'=getwd(),'nproc'=${task.cpus},'outname'='all_reps','log'='all_reps_clone_command_log'))"
    mv enchantr all_reps_dist_report
    """
}
