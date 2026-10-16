/*
 * Reformat design file and check validity
 */
process VALIDATE_INPUT {
    tag "$samplesheet"
    label 'immcantation'
    label 'process_single'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"

    input:
    file samplesheet
    path miairr
    val collapseby
    val cloneby
    val reassign

    output:
    path "*/validated_input.tsv", emit: validated_input
    path "*/validated_input_not-valid.tsv", emit: not_valid_input, optional: true
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    """
    Rscript -e "enchantr:::enchantr_report('validate_input', report_params=list('input'='${samplesheet}','collapseby'='${collapseby}','cloneby'='${cloneby}','reassign'='${reassign}','miairr'='${miairr}','outdir'=getwd()))"

    cp -r enchantr validate_input_report && rm -rf enchantr

    """
}
