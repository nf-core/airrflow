process FILTER_QUALITY {
    tag "$meta.id"
    label 'immcantation'
    label 'single_cpu'

    // TODO: update container
    container "immcantation/suite:devel"

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
