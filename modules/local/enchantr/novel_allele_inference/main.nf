def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process NOVEL_ALLELE_INFERENCE {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // meta, sequence tsv in AIRR format, reference fasta

    output:
    tuple val(meta), path("*_report/db_novel"), emit: reference // reference folder
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report", optional: true, emit: report
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('tigger'), eval("Rscript -e \"library(tigger); cat(as.character(packageVersion('tigger')))\""), emit: versions_tigger, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = tabs.join(',')
    """
    Rscript -e "enchantr::enchantr_report('novel_allele_inference', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus}, \\
                                        'log'='${meta.id}_novel_allele_inference_command_log' ${args}))"

    cp -r enchantr ${meta.id}_novel_allele_inference_report && rm -rf enchantr

    """
}
