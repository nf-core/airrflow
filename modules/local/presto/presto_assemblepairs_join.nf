process PRESTO_ASSEMBLEPAIRS_JOIN {
    tag "$meta.id"
    label 'process_long_parallelized'
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2), path(reads_pass)

    output:
    tuple val(meta), path("*_assemblejoin-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*.log")
    path("*_table.tab")
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    AssemblePairs.py join -1 $R1 -2 $R2 --nproc ${task.cpus} \\
        $args \\
        --outname ${meta.id}_join --log ${meta.id}_join.log > ${meta.id}_join_command_log.txt
    ParseLog.py -l ${meta.id}_join.log $args2
    cp ${meta.id}_assemble-pass.fastq ${meta.id}_assemblejoin-pass.fastq
    cat ${meta.id}_join_assemble-pass.fastq >> ${meta.id}_assemblejoin-pass.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( AssemblePairs.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
