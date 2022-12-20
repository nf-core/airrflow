process PRESTO_PARSEHEADERS_METADATA {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::presto=0.7.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ParseHeaders.py add -s $reads -o "${reads.baseName}_reheader-pass.fastq" $args -u ${meta.id} ${meta.subject_id} ${meta.species} ${meta.locus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( ParseHeaders.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
