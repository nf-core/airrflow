process PRESTO_MASKPRIMERS_SCORE {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.4--pyhdfd78af_0' :
        'biocontainers/presto:0.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(read)
    path(primers)
    val(primer_start)
    val(barcode)
    val(primer_maxerror)
    val(primer_mask_mode)
    val(reverse_primers)
    val(suffix)

    output:
    tuple val(meta), path("*_primers-pass.fastq"), emit: reads
    path "*.txt", emit: logs
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions


    script:
    def args = task.ext.args?: ''
    def args2 = task.ext.args2?: ''
    def revpr = reverse_primers ? '--revpr' : ''
    def barcode = barcode ? '--barcode' : ''
    """
    MaskPrimers.py score \\
    --nproc ${task.cpus} \\
    -s $read \\
    -p ${primers} \\
    --maxerror ${primer_maxerror} \\
    --mode ${primer_mask_mode} \\
    --start ${primer_start} \\
    $barcode $revpr \\
    $args \\
    --outname ${meta.id}_${suffix} \\
    --log ${meta.id}_${suffix}.log > ${meta.id}_command_log_${suffix}.txt
    ParseLog.py -l *.log $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
