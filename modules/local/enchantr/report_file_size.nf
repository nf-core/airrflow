/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'single_cpu'
    label 'enchantr'

    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.1--r41hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.0.1--r41hdfd78af_0' }"

    input:
    path logs

    output:
    path "*_report", emit: file_size
    path "versions.yml", emit: versions

    script:
    """
    echo "${logs.join('\n')}" > logs.txt
    Rscript -e "enchantr::enchantr_report('file_size', \\
        report_params=list('input'='logs.txt',\\
        'outdir'=getwd()))"

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr file_size_report
    """
}
