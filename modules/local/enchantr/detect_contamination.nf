process DETECT_CONTAMINATION {
    tag "multi_repertoire"

    label 'process_long_parallelized'
    label 'immcantation'


    conda "bioconda::r-enchantr=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.1--r42hdfd78af_0':
        'biocontainers/r-enchantr:0.1.1--r42hdfd78af_0' }"

    input:
    path(tabs)

    output:
    path("*cont-flag.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml" , emit: versions

    script:
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('contamination', \\
        report_params=list('input'='tabs.txt',\\
        'input_id'='id','outdir'=getwd(), \\
        'outname'='cont-flag', \\
        'log'='all_reps_contamination_command_log'))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    mv enchantr all_reps_cont_report
    """
}
