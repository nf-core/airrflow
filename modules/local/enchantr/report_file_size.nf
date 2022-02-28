/*
 * Generate file size report
 */
process REPORT_FILE_SIZE {
    tag "file_size"
    label 'immcantation'
    label 'enchantr'
    label 'single_cpu'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    path logs

    output:
    path "*_report", emit: file_size

    script:
    """
    echo "${logs.join('\n')}" > logs.txt
    Rscript -e "enchantr::enchantr_report('file_size', report_params=list('input'='logs.txt','outdir'=getwd()))"
    mv enchantr file_size_report
    """
}
