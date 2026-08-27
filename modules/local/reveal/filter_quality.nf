process FILTER_QUALITY {
    tag "$meta.id"
    label 'immcantation'
    label 'process_single'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*quality-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('alakazam'), eval('Rscript -e "library(alakazam); cat(paste(packageVersion(\'alakazam\'), collapse=\'.\'))"'), emit: versions_alakazam, topic: versions
    tuple val("${task.process}"), val('optparse'), eval('Rscript -e "library(optparse); cat(paste(packageVersion(\'optparse\'), collapse=\'.\'))"'), emit: versions_optparse, topic: versions
    tuple val("${task.process}"), val('stringi'), eval('Rscript -e "library(stringi); cat(paste(packageVersion(\'stringi\'), collapse=\'.\'))"'), emit: versions_stringi, topic: versions
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "library(dplyr); cat(paste(packageVersion(\'dplyr\'), collapse=\'.\'))"'), emit: versions_dplyr, topic: versions
    tuple val("${task.process}"), val('airr'), eval('Rscript -e "library(airr); cat(paste(packageVersion(\'airr\'), collapse=\'.\'))"'), emit: versions_airr, topic: versions
    tuple val("${task.process}"), val('DT'), eval('Rscript -e "library(DT); cat(paste(packageVersion(\'DT\'), collapse=\'.\'))"'), emit: versions_dt, topic: versions
    tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_r, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }

    """
    reveal_filter_quality.R --repertoire $tab --outname ${meta.id} > ${meta.id}_fq_command_log.txt

    """
}
