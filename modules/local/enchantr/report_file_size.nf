// Import generic module functions
include { saveFiles } from '../functions'

params.options = [:]

/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'enchantr'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report_file_size', publish_id:'') }

    conda (params.enable_conda ? { exit 1 "TODO: set up conda " } : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/python:3.8.3"  // Docker image
    }

    input:
    path logs

    output:
    path "enchantr", emit: file_size

    script:
    """
    Rscript -e "enchantr::enchantr_report('file_size', report_params=list('input'='${logs}','outdir'=getwd()))"
    """
}
