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

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:4.1.0':
        'docker.io/immcantation/airrflow:4.1.0' }"


    input:
    path tab // sequence tsv in AIRR format
    path logo
    path tabs_samplesheet

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "all_reps_dist_report/tables/*_threshold-summary.tsv", emit: threshold_summary, optional:true
    path "all_reps_dist_report/tables/*_threshold-mean.tsv", emit: mean_threshold
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ? asString(task.ext.args) : ''
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', \\
        report_params=list('input'='${tabs_samplesheet}',\\
            'cloneby'='${params.cloneby}',\\
            'crossby'='${params.crossby}',\\
            'singlecell'='${params.singlecell}',\\
            'outdir'=getwd(),\\
            'nproc'=${task.cpus},\\
            'outname'='all_reps',\\
            'log'='all_reps_threshold_command_log',\\
            'logo'='${logo}' ${args}))"

    cp -r enchantr all_reps_dist_report && rm -rf enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
