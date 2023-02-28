process PRESTO_COLLAPSESEQ {
    tag "$meta.id"
    label "process_medium"
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_collapse-unique.fastq") , emit: reads
    path("*_command_log.txt") , emit: logs
    path("*.log")
    path("*_table.tab")
    path("versions.yml"), emit: versions



    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    CollapseSeq.py -s $reads $args --outname ${meta.id} --log ${meta.id}.log > "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}.log" $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( CollapseSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
