include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FETCH_DATABASES {
    tag "IMGT IGBLAST"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::wget=1.20.1" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/wget:1.20.1"  // Singularity image
    } else {
        container "quay.io/biocontainers/wget:1.20.1"                        // Docker image
    }

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: imgt
    
    script:
    """
    echo "Fetching databases..."

    wget https://raw.githubusercontent.com/nf-core/test-datasets/bcellmagic/database-cache/databases.zip

    unzip databases.zip

    echo "FetchDBs process finished."
    """
}