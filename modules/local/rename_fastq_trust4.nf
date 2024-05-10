// Import generic module functions
process RENAME_FASTQ_TRUST4 {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.8.0 conda-forge::biopython=1.74"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' :
        'biocontainers/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' }"

    input:
    tuple val(meta), path(R1), path(R2)
    tuple val(meta_2), path(orig_r1), path(orig_r2)

    output:
    tuple val(meta), path(orig_r1), path(orig_r2) , emit: reads

    script:
    """
    mv ${R1} ${orig_r1}
    mv ${R2} ${orig_r2}
    """
}
