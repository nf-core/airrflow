def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param].toString()
        value = value.isNumber() ? value : "'${value}'"
        ",'${param}'=${value}"
    }.join('')
}

process COLLAPSE_DUPLICATES {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"

    input:
    tuple val(meta), path(tabs) // tuple [val(meta), compressed sequence tsv in AIRR format ]
    val collapseby

    output:
    tuple val(meta), path("*/*/*collapse-pass.tsv.gz"), emit: tab // compressed sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('collapse_duplicates', \\
        report_params=list('input'='tabs.txt',\\
        'collapseby'='${collapseby}',\\
        'outdir'=getwd(),\\
        'nproc'=${task.cpus},\\
        'outname'='${meta.id}',\\
        'log'='${meta.id}_collapse_command_log' ${args}))"

    cp -r enchantr ${meta.id}_collapse_report && rm -r enchantr

    """
}
