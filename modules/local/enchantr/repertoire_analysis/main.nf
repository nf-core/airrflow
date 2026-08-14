def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process REPERTOIRE_ANALYSIS {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs) // meta, sequence tsv in AIRR format
    path repertoires_samplesheet
    val cloneby

    output:
    tuple val(meta), path("*/*/*repertoire-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions


    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = ""
    if (repertoires_samplesheet) {
        input = repertoires_samplesheet
    } else {
        input = tabs.join(',')
    }
    """
    Rscript -e "enchantr::enchantr_report('repertoire_analysis', \\
                                        report_params=list('input'='${input}', \\
                                        'cloneby'='${cloneby}', \\
                                        'outputby'='${cloneby}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus}, \\
                                        'log'='${meta.id}_clone_command_log' ${args}))"

    cp -r enchantr repertoire_analysis_report && rm -rf enchantr

    """
}
