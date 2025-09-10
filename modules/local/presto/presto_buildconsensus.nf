process PRESTO_BUILDCONSENSUS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/presto:0.7.6--pyhdfd78af_0' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*_R1_consensus-pass.fastq"), path("*_R2_consensus-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*.tab"), emit: log_tab
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    """
    BuildConsensus.py -s $R1 --nproc ${task.cpus} ${args} --outname ${meta.id}_R1 --log ${meta.id}_R1.log > ${meta.id}_command_log.txt
    BuildConsensus.py -s $R2 --nproc ${task.cpus} ${args2} --outname ${meta.id}_R2 --log ${meta.id}_R2.log >> ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log ${args3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( BuildConsensus.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
