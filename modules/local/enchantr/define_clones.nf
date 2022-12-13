def asString (args) {
    s = ""
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

process DEFINE_CLONES {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.5--r42hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.0.5--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(tabs) // meta, sequence tsv in AIRR format
    val threshold
    path imgt_base

    output:
    path("*/*/*clone-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml", emit: versions


    script:
    def args = asString(task.ext.args) ?: ''
    def thr = threshold.join("")
    """
    Rscript -e "enchantr::enchantr_report('define_clones', \\
                                        report_params=list('input'='${tabs.join(',')}', \\
                                        'imgt_db'='${imgt_base}', \\
                                        'cloneby'='${params.cloneby}', \\
                                        'force'=FALSE, \\
                                        'threshold'=${thr}, \\
                                        'singlecell'='${params.singlecell}','outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='${meta.id}_clone_command_log' ${args}))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr '${meta.id}_clone_report'
    """
}
