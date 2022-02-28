process FILTER_JUNCTION_MOD3 {
    tag "$meta.id"
    label 'immcantation'
    label 'single_cpu'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*junction-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    """
    reveal_mod_3_junction.R --repertoire $tab --outname ${meta.id} > "${meta.id}_jmod3_command_log.txt"
    """
}
