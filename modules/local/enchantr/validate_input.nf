// Import generic module functions
include { saveFiles } from '../functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process VALIDATE_INPUT {
    tag "$samplesheet"
    label 'immcantation'
    label 'enchantr'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'validated_input', publish_id:'') }

    conda (params.enable_conda ? { exit 1 "TODO: set up conda for SAMPLESHEET_CHECK_AIRR." } : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/python:3.8.3"  // Docker image
    }

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
