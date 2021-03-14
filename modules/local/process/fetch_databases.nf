include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FETCH_DATABASES {
    tag "IMGT IGBLAST"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-91e61e7331dcf9dede9440679ae179eef1d60410:d7c34e324b30d1d37f98bb56af08b91a4725aa2a-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-91e61e7331dcf9dede9440679ae179eef1d60410:d7c34e324b30d1d37f98bb56af08b91a4725aa2a-0"                        // Docker image
    }

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: imgt
    
    script:
    //TODO: get db versions
    """
    echo "Fetching databases..."

    fetch_imgt.sh -o imgtdb_base

    fetch_igblastdb.sh -x -o igblast_base

    imgt2igblast.sh -i ./imgtdb_base -o igblast_base

    echo "FetchDBs process finished."
    """
}