process ADD_META_TO_TAB {
    tag "$meta.id"
    label 'immcantation'
    label 'process_single'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"

    cache 'deep' // Without 'deep' this process would run when using -resume

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(validated_input)

    output:
    tuple val(meta), path("*meta-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "library(dplyr); cat(paste(packageVersion(\'dplyr\'), collapse=\'.\'))"'), emit: versions_dplyr, topic: versions
    tuple val("${task.process}"), val('optparse'), eval('Rscript -e "library(optparse); cat(paste(packageVersion(\'optparse\'), collapse=\'.\'))"'), emit: versions_optparse, topic: versions
    tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_r, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    """
    reveal_add_metadata.R --repertoire ${tab} --metadata ${validated_input} --input_id ${meta.id} --outname ${meta.id} > ${meta.id}_add-meta_command_log.txt

    """
}
