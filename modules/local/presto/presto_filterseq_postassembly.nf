process PRESTO_FILTERSEQ_POSTASSEMBLY {
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

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
