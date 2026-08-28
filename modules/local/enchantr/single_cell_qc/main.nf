def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process SINGLE_CELL_QC {
    tag 'all_single_cell'
    label 'immcantation'
    label 'process_medium'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    path(tabs)

    output:
    path("*/*/*scqc-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path("*_report"), emit: report
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('single_cell_qc', \\
        report_params=list('input'='tabs.txt',\\
        'outdir'=getwd(), \\
        'log'='all_reps_scqc_command_log'  ${args} ))"

    cp -r enchantr all_reps_scqc_report && rm -rf enchantr

    """
}
