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
    tag "$tabs"
    label 'immcantation'
    label 'enchantr'
    label 'process_long'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    //tuple val(meta), path(tabs) // sequence tsv in AIRR format
    path(tabs)

    output:
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"

    script:
    meta=[]
    def args = asString(task.ext.args)
    def id_name = "$tabs".replaceFirst('__.*','')
    // TODO use nice outname, not tabs
    """
    Rscript -e "enchantr::enchantr_report('dowser_lineage', \\
                                        report_params=list('input'='${tabs}', \\
                                        'outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='${id_name}_dowser_command_log' ${args}))"
    mv enchantr '${id_name}_dowser_report'
    """
}
