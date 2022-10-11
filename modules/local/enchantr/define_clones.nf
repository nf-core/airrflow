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

process DEFINE_CLONES {
    tag 'all_reps'
    label 'immcantation'
    label 'enchantr'
    label 'process_high'
    label 'process_long'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    //tuple val(meta), path(tabs) // sequence tsv in AIRR format
    path(tabs)
    val threshold
    path imgt_base

    output:
    path("*clone-pass.tsv"), emit: tab, optional: true // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"

    script:
    def outname = ''
    if (task.ext.args.containsKey('outname')) { outname = task.ext.args['outname'] }
    """
    Rscript -e "enchantr::enchantr_report('define_clones', \\
                                        report_params=list('input'='${tabs.join(',')}', \\
                                        'imgt_db'='${imgt_base}', \\
                                        'cloneby'='${params.cloneby}','threshold'=${threshold}, \\
                                        'outputby'='id', \\
                                        'outname'='${outname}', \\
                                        'singlecell'='${singlecell}','outdir'=getwd(), \\
                                        'nproc'=${task.cpus},\\
                                        'log'='all_reps_clone_command_log' ${args}))"
    mv enchantr 'all_reps_clone_report'
    """
}
