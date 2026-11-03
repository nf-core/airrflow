process PRESTO_MASKPRIMERS_ALIGN {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1)
    path(primers)
    val(max_len)
    val(barcode)
    val(max_error)
    val(mask_mode)
    val(reverse_primers)
    val(suffix)

    output:
    tuple val(meta), path("*_primers-pass.fastq") , emit: reads
    path "*.txt", emit: logs
    path "*.tab", emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('MaskPrimers.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def barcode_param = barcode ? '--barcode' : ''
    def revpr = reverse_primers ? '--revpr' : ''
    def args = task.ext.args?: ''
    def args2 = task.ext.args2?: ''
    """
    MaskPrimers.py align --nproc ${task.cpus} \\
    -s $R1 \\
    -p ${primers} \\
    --maxlen ${max_len} \\
    --maxerror ${max_error} \\
    --mode ${mask_mode} \\
    $barcode_param \\
    $revpr \\
    $args \\
    --outname ${meta.id}_${suffix} \\
    --log ${meta.id}_${suffix}.log > ${meta.id}_command_log_${suffix}.txt
    ParseLog.py -l ${meta.id}_${suffix}.log $args2

    """
}
