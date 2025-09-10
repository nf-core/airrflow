process PRESTO_PARSEHEADERS {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/presto:0.7.6--pyhdfd78af_0' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads
    path "versions.yml" , emit: versions

    script:
    def subcommand = task.ext.subcommand?: ''
    def args = task.ext.args?: ''
    """
    ParseHeaders.py $subcommand -s $reads -o ${reads.baseName}_reheader-pass.fastq $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( ParseHeaders.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
