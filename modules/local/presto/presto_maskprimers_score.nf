process PRESTO_MASKPRIMERS_SCORE {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c538a0e310c303233164cbe486b5e5f6bddcf18975a9b20ac2f590f151f03e62/data' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

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
