

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_TABLES {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: tab // sequence tsv in AIRR format

    script:
    """
    echo "${meta.id}"
    echo "${meta.samples}"
    echo "${tab}"
    echo "${tab.join('\n')}" > tab.list
    
    head -n 1 ${tab[0]} > ${meta.id}.tsv
    tail -n +2 ${tab} >> ${meta.id}.tsv
    sed -i '/==>/d' ${meta.id}.tsv
    """
}
