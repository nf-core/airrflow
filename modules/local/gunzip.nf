process GUNZIP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("${R1.simpleName}*"), path("${R2.simpleName}*")   , emit: reads
    path "versions.yml", emit: versions

    script:
    """
    gunzip -f "${R1}"
    gunzip -f "${R2}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
