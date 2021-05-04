// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUNZIP {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path('*_R1.fastq'), path('*_R2.fastq')   , emit: reads
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    gunzip -f "${R1}"
    gunzip -f "${R2}"
    echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//' > ${software}.version.txt
    """
}