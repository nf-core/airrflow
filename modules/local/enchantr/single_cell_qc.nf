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

process SINGLE_CELL_QC {
    tag 'all_single_cell'
    label 'immcantation'
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:3.3.1':
        'docker.io/immcantation/airrflow:3.3.1' }"

    input:
    path(tabs)

    output:
    path("*/*/*scqc-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path("*_report"), emit: report
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ? asString(task.ext.args) : ''
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('single_cell_qc', \\
        report_params=list('input'='tabs.txt',\\
        'outdir'=getwd(), \\
        'log'='all_reps_scqc_command_log'  ${args} ))"

    cp -r enchantr all_reps_scqc_report && rm -rf enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
