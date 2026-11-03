process PRESTO_CLUSTERSETS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*_R1_cluster-pass.fastq"), path("*_R2_cluster-pass.fastq"), emit: reads
    path "*_command_log.txt", emit: logs
    path "*.tab", emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('ClusterSets.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version &> vsearch.txt; cat vsearch.txt | head -n 1 | grep -o \'v[0-9\\.]\\+\''), emit: versions_vsearch, topic: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    ClusterSets.py set --nproc ${task.cpus} -s $R1 --outname ${meta.id}_R1 $args --log ${meta.id}_R1.log > ${meta.id}_command_log.txt
    ClusterSets.py set --nproc ${task.cpus} -s $R2 --outname ${meta.id}_R2 $args --log ${meta.id}_R2.log >> ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log $args2

    """
}
