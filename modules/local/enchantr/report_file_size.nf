/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:devel':
        'docker.io/immcantation/airrflow:devel' }"

    input:
    path logs
    path metadata
    path logs_tabs

    output:
    path "*_report", emit: file_size
    path "versions.yml", emit: versions
    path "file_size_report/tables/log_data.tsv", emit: table

    script:
    """
    Rscript -e "enchantr::enchantr_report('file_size', \\
        report_params=list('input'='${logs_tabs}', 'metadata'='${metadata}',\\
        'outdir'=getwd()))"

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr file_size_report
    """
}
