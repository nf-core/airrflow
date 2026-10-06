process UNZIP_DB {
    tag "unzip_db"
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c8501d51f053cc82a0f87fe60516b362df73b37e56175d57d526f300684b802/data' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(archive)

    output:
    path("$unzipped")   , emit: unzipped
    tuple val("${task.process}"), val('unzip'), eval('unzip -v 2>&1 | head -n 1 | sed \'s/^.*UnZip //; s/ of.*$//\''), emit: versions_unzip, topic: versions

    script:
    unzipped = archive.toString() - '.zip'
    """
    unzip $archive

    """
}
