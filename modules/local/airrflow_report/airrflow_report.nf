process AIRRFLOW_REPORT {
    tag "${meta.id}"
    label 'process_high'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    path("Table_sequences.tsv")
    path("Table_sequences_assembled.tsv")
    path(repertoire_report)
    path(css)
    path(logo)
    path("reads_per_UMI.csv")

    output:
    tuple val("${task.process}"), val('alakazam'), eval("Rscript -e \"library(alakazam); cat(as.character(packageVersion('alakazam')))\""), emit: versions_alakazam, topic: versions
    tuple val("${task.process}"), val('shazam'), eval("Rscript -e \"library(shazam); cat(as.character(packageVersion('shazam')))\""), emit: versions_shazam, topic: versions
    tuple val("${task.process}"), val('stringr'), eval("Rscript -e \"library(stringr); cat(as.character(packageVersion('stringr')))\""), emit: versions_stringr, topic: versions
    tuple val("${task.process}"), val('dplyr'), eval("Rscript -e \"library(dplyr); cat(as.character(packageVersion('dplyr')))\""), emit: versions_dplyr, topic: versions
    tuple val("${task.process}"), val('knitr'), eval("Rscript -e \"library(knitr); cat(as.character(packageVersion('knitr')))\""), emit: versions_knitr, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e \"cat(as.character(getRversion()))\""), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('scales'), eval("Rscript -e \"library(scales); cat(as.character(packageVersion('scales')))\""), emit: versions_scales, topic: versions
    path("repertoire_comparison"), emit: results_folder
    path("*.html"), emit: report_html

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    """
    execute_report.R --report_file ${repertoire_report}

    """
}
