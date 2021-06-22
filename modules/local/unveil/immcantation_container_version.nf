// Import generic module functions
include { saveFiles; getSoftwareName } from '../functions'

params.options = [:]

/*
 * Immcantation version 
 */
process IMMCANTATION {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'immcantation_version', publish_id:'') }
    if (params.immcantation_container) {
        container params.immcantation_container
    }
    output:
    path "*.version.txt", emit: version
       
    script:
    def software = getSoftwareName(task.process)
    """
    versions report | head -n 1 > ${software}.version.txt
    """
}