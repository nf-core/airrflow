process UNZIP_CELLRANGERDB {
    tag "unzip_cellrangerdb"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(archive)

    output:
    path("$unzipped")   , emit: unzipped
    tuple val("${task.process}"), val('unzip_cellrangerdb'), eval('tar --version 2>&1 | head -n 1 | sed \'s/^.*(GNU tar) //; s/ Copyright.*$//\''), emit: versions_unzip_cellrangerdb, topic: versions

    script:
    unzipped = archive.toString() - '.tar.gz'
    """
    echo "${unzipped}"

    tar -xzvf ${archive}

    """
}
