def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process CLONAL_ASSIGNMENT {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // meta, sequence tsv in AIRR format
    val threshold
    path repertoires_samplesheet
    val cloneby
    val singlecell

    output:
    tuple val(meta), path("*/*/*clone-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('scoper'), eval("Rscript -e \"library(scoper); cat(as.character(packageVersion('scoper')))\""), emit: versions_scoper, topic: versions


    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def thr = threshold.join("")
    def input = ""
    if (repertoires_samplesheet) {
        input = repertoires_samplesheet
    } else {
        input = tabs.join(',')
    }
    """
    Rscript -e "enchantr::enchantr_report('clonal_assignment', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'species'='auto', \\
                                        'cloneby'='${cloneby}', \\
                                        'outputby'='${cloneby}', \\
                                        'force'=FALSE, \\
                                        'threshold'=${thr}, \\
                                        'singlecell'='${singlecell}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus}, \\
                                        'log'='${meta.id}_clone_command_log' ${args}))"

    cp -r enchantr ${meta.id}_clone_report && rm -rf enchantr

    """
}
