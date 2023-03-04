/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'single_cpu'

    conda "bioconda::r-enchantr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.0--r42hdfd78af_2':
        'quay.io/biocontainers/r-enchantr:0.1.0--r42hdfd78af_2' }"

    input:
    path logs
    path metadata

    output:
    path "*_report", emit: file_size
    path "versions.yml", emit: versions
    path "file_size_report/tables/log_data.tsv", emit: table

    script:
    """
    echo "${logs.join('\n')}" > logs.txt
    Rscript -e "enchantr::enchantr_report('file_size', \\
        report_params=list('input'='logs.txt', 'metadata'='${metadata}',\\
        'outdir'=getwd()))"

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr file_size_report
    """
}
