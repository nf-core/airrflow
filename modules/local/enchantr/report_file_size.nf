/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'process_single'

    conda "bioconda::r-enchantr=0.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/ssnn/suite:prerelease':
        'docker.io/ssnn/suite:prerelease' }"

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
