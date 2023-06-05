process CHANGEO_CONVERTDB_FASTA {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.0 bioconda::igblast=1.19.0 conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:00534555924705cdf2f7ac48b4b8b4083527ca58-1' :
        'biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:00534555924705cdf2f7ac48b4b8b4083527ca58-1' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*.fasta"), emit: fasta // sequence tsv in AIRR format
    path "versions.yml" , emit: versions
    path "*_command_log.txt" , emit: logs

    script:
    def args = task.ext.args ?: ''
    """
    ConvertDb.py fasta -d $tab $args > ${meta.id}_convertdb_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( ConvertDb.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
