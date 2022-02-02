// Import generic module functions
include { saveFiles; getSoftwareName } from '../functions'

params.options = [:]

/*
 * Immcantation version
 */
process IMMCANTATION {
    label 'immcantation'
    label 'single_cpu'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'immcantation_version', publish_id:'') }

    output:
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    if ! command -v versions report &> /dev/null
    then
        echo "immcantation: none" > ${software}.version.txt
    else
        versions report | head -n 1 > ${software}.version.txt
    fi
    """
}
