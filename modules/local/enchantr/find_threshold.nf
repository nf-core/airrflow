def asString (args) {
    def s = ""
    def value = ""
    if (args.size()>0) {
        if (args[0] != 'none') {
            for (param in args.keySet().sort()){
                value = args[param].toString()
                if (!value.isNumber()) {
                    value = "'"+value+"'"
                }
                s = s + ",'"+param+"'="+value
            }
        }
    }
    return s
}

process FIND_THRESHOLD {
    tag "all_reps"

    label 'process_long_parallelized'
    label 'immcantation'

    conda "bioconda::r-enchantr=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.1--r42hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.1.1--r42hdfd78af_0' }"


    input:
    path tab // sequence tsv in AIRR format
    path logo

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "*_threshold-summary.tsv", emit: threshold_summary
    path "*_threshold-mean.tsv", emit: mean_threshold
    path "versions.yml", emit: versions

    script:
    def args = asString(task.ext.args) ?: ''
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', \\
        report_params=list('input'='${tab.join(',')}',\\
            'cloneby'='${params.cloneby}',\\
            'crossby'='${params.crossby}',\\
            'singlecell'='${params.singlecell}',\\
            'outdir'=getwd(),\\
            'nproc'=${task.cpus},\\
            'outname'='all_reps',\\
            'log'='all_reps_threshold_command_log',\\
            'logo'='${logo}' ${args}))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    mv enchantr all_reps_dist_report
    """
}
