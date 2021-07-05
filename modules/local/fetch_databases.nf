include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FETCH_DATABASES {
    tag "IMGT IGBLAST"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'databases', publish_id:'') }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    if (!params.custom_container) {
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e5bd866a6803c301bf10b82e697c5dd2e49810a1-1"  // Singularity image
        } else {
            container "quay.io/biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e5bd866a6803c301bf10b82e697c5dd2e49810a1-1"                        // Docker image
        }
    } else {
        container params.custom_container
    }

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: imgt
    path "*.version.txt" , emit: version
    
    script:
    //TODO: get db versions. Currently using the download date
    def software = getSoftwareName(task.process)
    """
    fetch_databases.sh
    sed -n '2p' imgtdb_base/IMGT.yaml | awk -F' '  '{print \$2}' > ${software}-IMGT.version.txt
    """
}