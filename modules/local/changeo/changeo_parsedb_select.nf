process CHANGEO_PARSEDB_SELECT {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda (params.enable_conda ? "bioconda::changeo=1.2.0 bioconda::igblast=1.17.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*parse-select.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    if (meta.locus == 'IG'){
        """
        ParseDb.py select -d $tab $args --outname ${meta.id} > "${meta.id}_select_command_log.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
            changeo: \$( ParseDb.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else if (meta.locus == 'TR'){
        """
        ParseDb.py select -d $tab $args2 --outname ${meta.id} > "${meta.id}_command_log.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
            changeo: \$( ParseDb.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    }
}
