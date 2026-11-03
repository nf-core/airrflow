process CHANGEO_CREATEGERMLINES {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    path(reference_fasta) // reference fasta

    output:
    tuple val(meta), path("*germ-pass.tsv"), emit: tab
    path("*_command_log.txt"), emit: logs
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('changeo'), eval('CreateGermlines.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    CreateGermlines.py -d ${tab} \\
    -r ${reference_fasta}/${meta.species}/vdj/ \\
    -g dmask --format airr \\
    --log ${meta.id}.log --outname ${meta.id} $args > ${meta.id}_create-germlines_command_log.txt
    ParseLog.py -l ${meta.id}.log -f ID V_CALL D_CALL J_CALL

    """
}
