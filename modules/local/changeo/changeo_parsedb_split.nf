process CHANGEO_PARSEDB_SPLIT {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/changeo_igblast_wget:dcfe290eb28df215' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*productive-T.tsv"), emit: tab //sequence tsv in AIRR format
    tuple val(meta), path("*productive-F.tsv"), emit: unproductive, optional: true //optional non-productive sequences
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
