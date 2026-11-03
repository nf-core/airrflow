process PRESTO_COLLAPSESEQ {
    tag "$meta.id"
    label "process_medium"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_collapse-unique.fastq") , emit: reads
    path("*_command_log.txt") , emit: logs
    path("*_table.tab")
    tuple val("${task.process}"), val('presto'), eval('CollapseSeq.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions



    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    CollapseSeq.py -s $reads $args --outname ${meta.id} --log ${meta.id}.log > ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}.log $args2

    """
}
