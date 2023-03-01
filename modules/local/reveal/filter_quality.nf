process FILTER_QUALITY {
    tag "$meta.id"
    label 'immcantation'
    label 'single_cpu'

    conda "bioconda::r-enchantr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.0--r42hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.1.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*quality-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    // TODO: add to enchantr
    """
    reveal_filter_quality.R --repertoire $tab --outname ${meta.id} > "${meta.id}_fq_command_log.txt"
    """
}
