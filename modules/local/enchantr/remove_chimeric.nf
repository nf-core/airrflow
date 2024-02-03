process REMOVE_CHIMERIC {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'


    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:3.3.0':
        'docker.io/immcantation/airrflow:3.3.0' }"


    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(imgt_base)

    output:
    tuple val(meta), path("*chimera-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" //, emit: chimera_report
    path "versions.yml" , emit: versions

    script:
    """
    Rscript -e "enchantr:::enchantr_report('chimera_analysis', \\
        report_params=list('input'='${tab}',\\
        'outdir'=getwd(), \\
        'nproc'=${task.cpus},\\
        'outname'='${meta.id}', \\
        'log'='${meta.id}_chimeric_command_log'))"

    cp -r enchantr ${meta.id}_chimera_report && rm -rf enchantr

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
