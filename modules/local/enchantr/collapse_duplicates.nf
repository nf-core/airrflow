process COLLAPSE_DUPLICATES {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:3.3.1':
        'docker.io/immcantation/airrflow:3.3.1' }"

    input:
    tuple val(meta), path(tabs) // tuple [val(meta), sequence tsv in AIRR format ]

    output:
    tuple val(meta), path("*/*/*collapse-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml" , emit: versions

    script:
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('collapse_duplicates', \\
        report_params=list('input'='tabs.txt',\\
        'collapseby'='${params.collapseby}',\\
        'outdir'=getwd(),\\
        'nproc'=${task.cpus},\\
        'outname'='${meta.id}',\\
        'log'='${meta.id}_collapse_command_log'))"

    cp -r enchantr ${meta.id}_collapse_report && rm -r enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
