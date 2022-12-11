def asString (args) {
    s = ""
    if (args.size()>0) {
        if (args[0] != 'none') {
            for (param in args.keySet().sort()){
                s = s + ",'"+param+"'='"+args[param]+"'"
            }
        }
    }
    return s
}

process DOWSER_LINEAGES {
    tag "${meta.id}"

    label 'process_high'
    label 'process_long'
    label 'error_ignore'
    label 'immcantation'
    label 'enchantr'

    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.5--r42hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.0.5--r42hdfd78af_0' }"

    input:
    //tuple val(meta), path(tabs) // sequence tsv in AIRR format
    tuple val(meta), path(tabs)

    output:
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml", emit: versions

    script:
    def args = asString(task.ext.args) ?: ''
    def id_name = "$tabs".replaceFirst('__.*','')
    // TODO use nice outname, not tabs
    """
    Rscript -e "enchantr::enchantr_report('dowser_lineage', \\
                                        report_params=list('input'='${tabs}', \\
                                        'exec'='${params.igphyml}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='${id_name}_dowser_command_log' ${args}))"

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr '${id_name}_dowser_report'
    """
}
