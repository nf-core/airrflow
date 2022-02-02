include { initOptions; saveFiles; getSoftwareName; asString } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process DOWSER_LINEAGES {
    tag "$tabs"
    label 'immcantation'
    label 'enchantr'
    label 'process_long'

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
    //tuple val(meta), path(tabs) // sequence tsv in AIRR format
    path(tabs)

    output:
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"

    script:
    meta=[]
    def args = asString(options.args)
    def id_name = "$tabs".replaceFirst('__.*','')
    // TODO use nice outname, not tabs
    """
    Rscript -e "enchantr::enchantr_report('dowser_lineage', \\
                                        report_params=list('input'='${tabs}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='${id_name}_dowser_command_log' ${args}))"
    mv enchantr '${id_name}_dowser_report'
    """
}
