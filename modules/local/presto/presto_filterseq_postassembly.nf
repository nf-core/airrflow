process PRESTO_FILTERSEQ_POSTASSEMBLY {
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
    tuple val(meta), path("*quality-pass.fastq") ,  emit: reads
    path "*_command_log.txt" , emit: logs
    path "versions.yml" , emit: versions
    path "*.log"
    path "*.tab" , emit: log_tab

    script:
    """
    FilterSeq.py quality -s $reads -q ${params.filterseq_q} --outname "${meta.id}" --log "${reads.baseName}.log" --nproc ${task.cpus} > "${meta.id}_command_log.txt"
    ParseLog.py -l "${reads.baseName}.log" -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
