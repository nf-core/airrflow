process PRESTO_PARSEHEADERS {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads

    script:
    def subcommand = task.ext.subcommand?: ''
    def args = task.ext.args?: ''
    """
    ParseHeaders.py $subcommand -s $reads -o "${reads.baseName}_reheader-pass.fastq" $args
    """
}
