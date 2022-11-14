process FIND_THRESHOLD {
    tag "all_reps"

    label 'process_high'
    label 'process_long'
    label 'enchantr'


    cache 'lenient'

    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.3--r42hdfd78af_1':
        'quay.io/biocontainers/r-enchantr:0.0.3--r42hdfd78af_1' }"


    input:
    path tab // sequence tsv in AIRR format

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "*_threshold-summary.tsv", emit: threshold_summary
    path "*_threshold-mean.tsv", emit: mean_threshold
    path "versions.yml", emit: versions

    script:
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', \\
        report_params=list('input'='${tab.join(',')}',\\
            'cloneby'='${params.cloneby}',\\
            'singlecell'='${params.singlecell}',\\
            'outdir'=getwd(),\\
            'nproc'=${task.cpus},\\
            'outname'='all_reps',\\
            'log'='all_reps_clone_command_log'))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr all_reps_dist_report
    """
}
