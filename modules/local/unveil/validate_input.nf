// Import generic module functions
include { saveFiles } from '../functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process VALIDATE_INPUT {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'validated_input', publish_id:'') }

    conda (params.enable_conda ? { exit 1 "TODO: set up conda for SAMPLESHEET_CHECK_AIRR." } : null)
    if (!params.custom_container) {
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
        } else {
            container "quay.io/biocontainers/python:3.8.3"  // Docker image
        }
    } else {
        container params.custom_container
    }

    input:
    file samplesheet
    path miairr
    val collapseby
    val cloneby
    
    output:
    path "validated_input.tsv", emit: validated_input
    path "validated_input_not-valid.tsv", emit: not_valid_input, optional: true
    path "validated_input.html", emit: validated_input_html

    script:  // This script is bundled with the pipeline, in bin/
    """
    unveil_validate_input.R --input "${samplesheet}" --collapseby ${collapseby} --cloneby ${cloneby} --output "validated_input" --miairr "${miairr}"
    """
}
