process CHANGEO_PARSEDB_SPLIT {
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
    tuple val(meta), path("*productive-T.tsv"), emit: tab //sequence tsv in AIRR format
    tuple val(meta), path("*productive-F.tsv"), emit: unproductive, optional: true //optional non-productive sequences
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('changeo'), eval('ParseDb.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions

    script:
    """
    ParseDb.py split -d $tab -f productive --outname ${meta.id} > ${meta.id}_split_command_log.txt

    """
}
