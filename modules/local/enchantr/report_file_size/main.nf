/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'process_single'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"

    input:
    path logs
    path metadata
    path logs_tabs

    output:
    path "*_report", emit: file_size
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    path "file_size_report/tables/log_data.tsv", emit: table

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    """
    Rscript -e "enchantr::enchantr_report('file_size', \\
        report_params=list('input'='${logs_tabs}', 'metadata'='${metadata}',\\
        'outdir'=getwd()))"

    cp -r enchantr file_size_report && rm -rf enchantr


    """
}
