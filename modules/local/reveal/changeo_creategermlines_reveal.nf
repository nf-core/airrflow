process CHANGEO_CREATEGERMLINES_REVEAL {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(imgt_base) // imgt db

    output:
    tuple val(meta), path("*germ-pass.tsv"), emit: tab
    path("*_command_log.txt"), emit: logs
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    CreateGermlines.py -d ${tab} -g dmask \\
    -r ${imgt_base}/${meta.species}/vdj/ --format airr --outdir . \\
    --log ${meta.id}.log --outname ${meta.id} $args > "${meta.id}_create-germlines_command_log.txt"
    ParseLog.py -l ${meta.id}.log -f ID V_CALL D_CALL J_CALL

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$( CreateGermlines.py --version | awk -F' '  '{print \$2}' )
        presto: \$( ParseLog.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}



