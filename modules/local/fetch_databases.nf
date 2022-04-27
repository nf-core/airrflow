process FETCH_DATABASES {
    tag "IMGT IGBLAST"
    label 'process_low'
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e5bd866a6803c301bf10b82e697c5dd2e49810a1-1' :
        'quay.io/biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e5bd866a6803c301bf10b82e697c5dd2e49810a1-1' }"

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: imgt
    path "versions.yml" , emit: versions

    script:
    """
    fetch_databases.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMGT download date: \$( echo \$(date "+%F") )
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
