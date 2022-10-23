process DETECT_CONTAMINATION {
    tag "multi_repertoire"

    label 'process_high'
    label 'process_long_parallelized'
    label 'immcantation'


    //conda (params.enable_conda ? "bioconda::r-enchantr=0.0.1" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.1--r41hdfd78af_0':
    //    'quay.io/biocontainers/r-enchantr:0.0.1--r41hdfd78af_0' }"
    container 'immcantation/suite:devel'
    // TODO: fix issue in enchantr missing r-reshape2 dependency and update container

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

    mv enchantr al_reps_cont_report
    """
}
