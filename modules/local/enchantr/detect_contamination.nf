process DETECT_CONTAMINATION {
    tag "multi_repertoire"

    label 'process_long_parallelized'
    label 'immcantation'


    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "docker.io/immcantation/airrflow:4.1.0"

    input:
    path(tabs)

    output:
    path("*cont-flag.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml" , emit: versions

    script:
    """
    echo "${tabs.join('\n')}" > tabs.txt
    Rscript -e "enchantr::enchantr_report('contamination', \\
        report_params=list('input'='tabs.txt',\\
        'input_id'='id','outdir'=getwd(), \\
        'outname'='cont-flag', \\
        'log'='all_reps_contamination_command_log'))"

    cp -r enchantr all_reps_cont_report && rm -rf enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
