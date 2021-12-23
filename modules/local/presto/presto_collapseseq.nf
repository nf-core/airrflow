process PRESTO_COLLAPSESEQ {
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_collapse-unique.fastq") , emit: reads
    path("*_command_log.txt") , emit: logs
    path("*.log")
    path("*_table.tab")


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    CollapseSeq.py -s $reads $args --outname ${meta.id} --log ${meta.id}.log > "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}.log" $args2
    """
}
