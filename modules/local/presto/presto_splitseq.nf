process PRESTO_SPLITSEQ {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::presto=0.7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.0--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_atleast-2.fasta"), emit: fasta
    path("*_command_log.txt"), emit: logs
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    SplitSeq.py group -s $reads \\
    $args \\
    --outname ${meta.id} \\
    --fasta > "${meta.id}_command_log.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( SplitSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
