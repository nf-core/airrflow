process PRESTO_SPLITSEQ {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/presto:0.7.6--ac08dbe217c927bd' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_atleast-*.fasta"), emit: fasta
    path("*_command_log.txt"), emit: logs
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    SplitSeq.py group -s $reads \\
    $args \\
    --outname ${meta.id} \\
    --fasta > ${meta.id}_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( SplitSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
