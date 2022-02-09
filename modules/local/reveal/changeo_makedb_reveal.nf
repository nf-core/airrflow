process CHANGEO_MAKEDB_REVEAL {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect()
    path(imgt_base)

    output:
    tuple val(meta), path("*db-pass.tsv"), emit: tab //sequence table in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    """
    MakeDb.py igblast -i $igblast -s $reads -r \\
    ${imgt_base}/${meta.species}/vdj/ \\
    $options.args \\
    --outname "${meta.id}" > "${meta.id}_mdb_command_log.txt"
    """
}
