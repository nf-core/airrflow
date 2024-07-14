process CHANGEO_PARSEDB_SPLIT {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.0 bioconda::igblast=1.22.0 conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' :
        'biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*productive-T.tsv"), emit: tab // sequence tsv in AIRR format
    tuple val(meta), path("*productive-F.tsv"), emit: unproductive
    path("*_command_log.txt"), emit: logs //process logs
    path "versions.yml" , emit: versions

    script:
    """
    ParseDb.py split -d $tab -f productive --outname ${meta.id} > ${meta.id}_split_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        changeo: \$( ParseDb.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
