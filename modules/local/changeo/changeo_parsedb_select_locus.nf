process CHANGEO_PARSEDB_SELECT_LOCUS {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.5 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0 bioconda::airr=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8d583137003a25ebc5fdd51c08fccacc3f3673c9fc162c92b006335a2cd6df0/data' :
        'community.wave.seqera.io/library/airr_changeo_igblast_wget:7303fb2eff967b2e' }"

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
