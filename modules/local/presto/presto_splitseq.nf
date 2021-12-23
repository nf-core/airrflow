process PRESTO_SPLITSEQ {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_atleast-2.fasta"), emit: fasta
    path("*_command_log.txt"), emit: logs

    script:
    """
    SplitSeq.py group -s $reads \\
    $options.args \\
    --outname ${meta.id} \\
    --fasta > "${meta.id}_command_log.txt"
    """
}
