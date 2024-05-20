// Import generic module functions
process RENAME_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.8.0 conda-forge::biopython=1.74"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' :
        'biocontainers/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}_${file.name}")  , emit: file

    script:
    """
    mv ${file} ${meta.id}.${file.extension}
    """
}
