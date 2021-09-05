include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process FILTER_QUALITY {
    tag "$meta.id"
    label 'immcantation'

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
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*quality-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    """
    reveal_filter_quality.R --repertoire $tab --outname ${meta.id} > "${meta.id}_${task.process}_command_log.txt"
    """
}
