process COLLAPSE_DUPLICATES {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'

    conda "bioconda::r-enchantr=0.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.3--r42hdfd78af_0':
        'biocontainers/r-enchantr:0.1.3--r42hdfd78af_0' }"

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

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr ${meta.id}_collapse_report
    """
}
