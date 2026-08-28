def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process DOWSER_LINEAGES {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs)
    val lineage_tree_builder
    val lineage_tree_exec

    output:
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('dowser'), eval("Rscript -e \"library(dowser); cat(as.character(packageVersion('dowser')))\""), emit: versions_dowser, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def id_name = "$tabs".replaceFirst('__.*','')
    // TODO use nice outname, not tabs
    """
    Rscript -e "enchantr::enchantr_report('dowser_lineage', \\
                                        report_params=list('input'='${tabs}', \\
                                        'build'='${lineage_tree_builder}', \\
                                        'exec'='${lineage_tree_exec}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='${id_name}_dowser_command_log' ${args}))"

    cp -r enchantr ${id_name}_dowser_report && rm -rf enchantr

    """
}
