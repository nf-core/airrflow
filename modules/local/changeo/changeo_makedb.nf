process CHANGEO_MAKEDB {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::changeo=1.2.0 bioconda::igblast=1.17.1" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' }"


    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect()
    path(imgt_base)

    output:
    tuple val(meta), path("*db-pass.tsv"), emit: tab //sequence table in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def args = task.ext.args ?: ''
    """
    MakeDb.py igblast -i $igblast -s $reads -r \\
    ${imgt_base}/${meta.species.toLowerCase()}/vdj/ \\
    $args \\
    --outname "${meta.id}" > "${meta.id}_command_log.txt"
    """
}
