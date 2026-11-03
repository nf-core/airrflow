process CHANGEO_PARSEDB_SELECT_LOCUS {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*parse-select.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('changeo'), eval('ParseDb.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions

    script:
    if (meta.locus.toUpperCase() == 'IG'){
        """
        ParseDb.py select -d $tab -f locus -u "IG[HKL]" --regex --outname ${meta.id} > ${meta.id}_select_command_log.txt

        """
    } else if (meta.locus.toUpperCase() == 'TR'){
        """
        ParseDb.py select -d $tab -f locus -u "TR[ABDG]" --regex --outname ${meta.id} > "${meta.id}_command_log.txt"

        """
    }
}
