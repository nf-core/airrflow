/*
 * Reformat design file and check validity
 */
process VALIDATE_INPUT {
    tag "$samplesheet"
    label 'immcantation'
    label 'enchantr'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    file samplesheet
    path miairr
    val collapseby
    val cloneby
    val reassign

    output:
    path "validated_input.tsv", emit: validated_input
    path "validated_input_not-valid.tsv", emit: not_valid_input, optional: true

    script:
    """
    Rscript -e "enchantr:::enchantr_report('validate_input', report_params=list('input'='${samplesheet}','collapseby'='${collapseby}','cloneby'='${cloneby}','reassign'='${reassign}','miairr'='${miairr}','outdir'=getwd()))"
    """
}
