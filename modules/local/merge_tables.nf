process MERGE_TABLES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

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
