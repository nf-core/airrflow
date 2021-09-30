// Import generic module functions
include { saveFiles; initOptions } from './functions'

params.options = [:]
options = initOptions(params.options)


process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    path samplesheet

    output:
    path '*.tsv'


    script:  // This script is bundled with the pipeline, in nf-core/dsltwotest/bin/
    """
    check_samplesheet.py $samplesheet
    cp $samplesheet samplesheet.valid.tsv
    """
}
