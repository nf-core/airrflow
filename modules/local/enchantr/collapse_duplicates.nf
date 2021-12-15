include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process COLLAPSE_DUPLICATES {
    tag "all_reps"
    label 'immcantation'
    label 'enchantr'
    label 'process_long'

    cache  'lenient'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
    }

    input:
    path(tabs) // tuple val(meta) // sequence tsv in AIRR format
    val(collapseby)

    output:
    tuple val(meta), path("*collapse-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report", emit: duplicates_report

    script:
    meta=[]
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('collapse_duplicates', report_params=list('input'='tabs.txt','collapseby'='${collapseby}','outdir'=getwd(), 'nproc'=${task.cpus},'outname'='all_reps', 'log'='all_reps_collapse_command_log'))"
    mv enchantr all_reps_collapse_report
    """
}
