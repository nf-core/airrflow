process PRESTO_MASKPRIMERS_EXTRACT {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(read)
    val(extract_start)
    val(extract_length)
    val(extract_mode)
    val(barcode)
    val(suffix)

    output:
    tuple val(meta), path("*_primers-pass.fastq") , emit: reads
    path "*.txt", emit: logs
    path "*.tab", emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('MaskPrimers.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args?: ''
    def args2 = task.ext.args2?: ''
    def barcode_param = barcode ? '--barcode' : ''
    """
    MaskPrimers.py extract \\
    --nproc ${task.cpus} \\
    -s $read \\
    --start ${extract_start} \\
    --len ${extract_length} \\
    $barcode_param \\
    $args \\
    --mode ${extract_mode} \\
    --outname ${meta.id}_${suffix} \\
    --log ${meta.id}_${suffix}.log >> ${meta.id}_command_log_${suffix}.txt
    ParseLog.py -l *.log $args2

    """
}
