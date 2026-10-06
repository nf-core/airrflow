process GUNZIP {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container 'docker.io/biocontainers/biocontainers:v1.2.0_cv1'

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("${R1.simpleName}*"), path("${R2.simpleName}*")   , emit: reads
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version 2>&1 | head -n 1 | sed \'s/^.*(gzip) //; s/ Copyright.*$//\''), emit: versions_gunzip, topic: versions

    script:
    """
    gunzip -f "${R1}"
    gunzip -f "${R2}"

    """
}
