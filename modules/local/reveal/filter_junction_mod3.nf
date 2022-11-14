process FILTER_JUNCTION_MOD3 {
    tag "$meta.id"
    label 'immcantation'
    label 'single_cpu'

    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.3--r42hdfd78af_1':
        'quay.io/biocontainers/r-enchantr:0.0.3--r42hdfd78af_1' }"

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
