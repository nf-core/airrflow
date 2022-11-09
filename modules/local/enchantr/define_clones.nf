process DEFINE_CLONES {
    tag 'all_reps'

    label 'process_high'
    label 'enchantr'


    //conda (params.enable_conda ? "bioconda::r-enchantr=0.0.1" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.1--r41hdfd78af_0':
    //    'quay.io/biocontainers/r-enchantr:0.0.1--r41hdfd78af_0' }"
    container 'immcantation/suite:devel'
    // TODO: fix issue in enchantr missing r-reshape2 dependency and update container

    input:
    //tuple val(meta), path(tabs) // sequence tsv in AIRR format
    path(tabs)
    val threshold
    path imgt_base

    output:
    path("*clone-pass.tsv"), emit: tab, optional: true // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml", emit: versions


    script:
    def outname = ''
    if (task.ext.args.containsKey('outname')) { outname = task.ext.args['outname'] }
    """
    Rscript -e "enchantr::enchantr_report('define_clones', \\
                                        report_params=list('input'='${tabs.join(',')}', \\
                                        'imgt_db'='${imgt_base}', \\
                                        'cloneby'='${params.cloneby}','threshold'=${threshold}, \\
                                        'outputby'='sample_id', \\
                                        'outname'='${outname}', \\
                                        'singlecell'='${params.singlecell}','outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='all_reps_clone_command_log' ${args}))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr 'all_reps_clone_report'
    """
}
