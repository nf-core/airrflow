process PRESTO_MASKPRIMERS {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(read)
    path(primers)
    val(barcode)
    val(primer_maxerror)
    val(primer_mask_mode)
    val(reverse_primers)
    val(suffix)

    output:
    tuple val(meta), path("*_primers-pass.fastq"), emit: reads
    path "*_command_log.txt", emit: logs
    path "*.log"
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions


    script:
    def revpr = params.reverse_primers ? '--revpr' : ''
    """
    MaskPrimers.py score --nproc ${task.cpus} -s $read -p ${primers} $barcode $revpr --maxerror ${primer_maxerror} --mode ${primer_mask_mode} --outname ${meta.id}_${suffix} --log ${meta.id}_${suffix}.log > ${meta.id}_${suffix}_command_log.txt
    ParseLog.py -l *.log -f ID PRIMER ERROR

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
