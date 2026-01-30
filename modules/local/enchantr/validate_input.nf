/*
 * Reformat design file and check validity
 */
process VALIDATE_INPUT {
    tag "$samplesheet"
    label 'immcantation'
    label 'process_single'
    label 'immcantation_container'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "docker.io/immcantation/airrflow:4.5.0dev"

    input:
    file samplesheet
    path miairr
    val collapseby
    val cloneby

    output:
    path "*/validated_input.tsv", emit: validated_input
    path "*/validated_input_not-valid.tsv", emit: not_valid_input, optional: true
    path "versions.yml", emit: versions

    script:
    """
    Rscript -e "enchantr:::enchantr_report('validate_input', report_params=list('input'='${samplesheet}','collapseby'='${collapseby}','cloneby'='${cloneby}','reassign'='${params.reassign}','miairr'='${miairr}','outdir'=getwd()))"

    cp -r enchantr validate_input_report && rm -rf enchantr

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
