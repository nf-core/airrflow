process PRESTO_MASKPRIMERS_ALIGN {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1)
    path(primers)
    val(max_len)
    val(max_error)
    val(mask_mode)
    val(suffix)

    output:
    tuple val(meta), path("*_primers-pass.fastq") , emit: reads
    path "*.txt", emit: logs
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args?: ''
    def args2 = task.ext.args2?: ''
    """
    MaskPrimers.py align --nproc ${task.cpus} \\
    -s $R1 \\
    -p ${primers} \\
    --maxlen ${max_len} \\
    --maxerror ${max_error} \\
    --mode ${mask_mode} \\
    $args \\
    --outname ${meta.id}_${suffix} \\
    --log ${meta.id}_${suffix}.log > ${meta.id}_command_log_${suffix}.txt
    ParseLog.py -l ${meta.id}_${suffix}.log $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
