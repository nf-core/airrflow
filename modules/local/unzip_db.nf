process UNZIP_DB {
    tag "unzip_db"
    label 'process_medium'

    conda "conda-forge::unzip=6.0 conda-forge::gzip=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4ba84b850752cf5e2d1d4857a30598f2939b01c59e25702b4c31ee547805490/data' :
        'community.wave.seqera.io/library/gzip_unzip:c40ea0e78704cb64' }"

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
