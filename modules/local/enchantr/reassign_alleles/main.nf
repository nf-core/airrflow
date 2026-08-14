def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process REASSIGN_ALLELES {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // meta, sequence tsv in AIRR format, reference fasta
    val segments // which segments to reassign alleles to
    val outputby // which field to use for output

    output:
    tuple val(meta), path("*/*/*reassign-pass.tsv"), emit: tab // reassigned repertoire
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions


    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def segs = segments.join(",")
    def input = tabs.join(',')

    """
    Rscript -e "enchantr::enchantr_report('reassign_alleles', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'species'='auto', \\
                                        'outputby'='${outputby}', \\
                                        'segments'='${segs}', \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_reassign_alleles_command_log' ${args}))"

    cp -r enchantr ${meta.id}_reassign_alleles_report && rm -rf enchantr

    """
}
