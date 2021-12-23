process CHANGEO_CONVERTDB_FASTA {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*.fasta"), emit: fasta // sequence tsv in AIRR format
    path "versions.yml" , emit: versions
    path "*_command_log.txt" , emit: logs

    script:
    def args = task.ext.args ?: ''
    """
    ConvertDb.py fasta -d $tab $args > "${meta.id}_command_log.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$( ConvertDb.py --version | awk -F' ' '{print \$2}' )
    END_VERSIONS
    """
}
