include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process FIND_THRESHOLD {
    tag "all_reps"
    label 'immcantation'
    label 'enchantr'
    label 'process_long'

    cache 'lenient'

    // TODO: update container
    container "immcantation/suite:devel"


    input:
    path tab // tuple val(meta) // sequence tsv in AIRR format
    val(cloneby)
    val(singlecell)

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "*_threshold-summary.tsv"
    path "*_threshold-mean.tsv", emit: mean_threshold

    script:
    meta=[]
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', report_params=list('input'='${tab.join(',')}','cloneby'='${cloneby}','singlecell'='${singlecell}','outdir'=getwd(), 'nproc'=${task.cpus},'outname'='all_reps', 'log'='all_reps_clone_command_log'))"
    mv enchantr all_reps_dist_report
    """
}
