// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path '*.tsv'


    script:  // This script is bundled with the pipeline, in nf-core/dsltwotest/bin/
    // TODO: update to check samplesheet
    """
    cp $samplesheet samplesheet.valid.tsv
    """
}