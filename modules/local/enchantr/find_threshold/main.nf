def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process FIND_THRESHOLD {
    tag "all_reps"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"


    input:
    path tab // sequence tsv in AIRR format
    path logo
    path tabs_samplesheet
    val cloneby
    val crossby
    val singlecell

    output:
    // tuple val(meta), path("*threshold-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "all_reps_dist_report/tables/*_threshold-summary.tsv", emit: threshold_summary, optional:true
    path "all_reps_dist_report/tables/*_threshold-mean.tsv", emit: mean_threshold
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions

    script:
        // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    """
    Rscript -e "enchantr::enchantr_report('find_threshold', \\
        report_params=list('input'='${tabs_samplesheet}',\\
            'cloneby'='${cloneby}',\\
            'crossby'='${crossby}',\\
            'singlecell'='${singlecell}',\\
            'outdir'=getwd(),\\
            'nproc'=${task.cpus},\\
            'outname'='all_reps',\\
            'log'='all_reps_threshold_command_log',\\
            'logo'='${logo}' ${args}))"

    cp -r enchantr all_reps_dist_report && rm -rf enchantr

    """
}
