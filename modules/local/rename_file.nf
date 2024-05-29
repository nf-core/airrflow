// Import generic module functions
process RENAME_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}.${file.extension}")  , emit: file

    script:
    """
    mv ${file} ${meta.id}.${file.extension}
    """
}
