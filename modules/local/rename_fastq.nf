// Import generic module functions
process RENAME_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.0 conda-forge::biopython=1.74" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' :
        'quay.io/biocontainers/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path('*_R1.fastq.gz'), path('*_R2.fastq.gz')   , emit: reads

    script:
    """
    mv "${R1}" "${meta.id}_R1.fastq.gz"
    mv "${R2}" "${meta.id}_R2.fastq.gz"
    """
}
