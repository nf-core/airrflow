process PRESTO_MASKPRIMERS_EXTRACT {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R2)
    val(extract_start),
    val(extract_length),
    val(extract_mode)

    output:
    tuple val(meta), path("*_R2_primers-pass.fastq") , emit: reads
    path "*_command_log_R2.txt", emit: logs
    path "*_R2.log"
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args?: ''
    def args2 = task.ext.args2?: ''
    """
    MaskPrimers.py extract --nproc ${task.cpus} \\
    -s $R2 \\
    --start ${extract_start} \\
    --len ${extract_length} \\
    $args \\
    --mode ${extract_mode} \\
    --outname ${meta.id}_R2 \\
    --log ${meta.id}_R2.log >> ${meta.id}_command_log_R2.txt
    ParseLog.py -l ${meta.id}_R2.log $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
