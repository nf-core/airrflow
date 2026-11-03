process CHANGEO_MAKEDB {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'


    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect()
    path(reference_fasta)

    output:
    tuple val(meta), path("*db-pass.tsv"), emit: tab //sequence table in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('changeo'), eval('MakeDb.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions

    script:
    def args = task.ext.args ?: ''
    def partial = meta.species.toLowerCase()=='mouse' && meta.locus.toLowerCase()=='tr'  ? '--partial' : ''
    """
    MakeDb.py igblast -i $igblast -s $reads -r \\
    ${reference_fasta}/${meta.species.toLowerCase()}/vdj/ \\
    --nproc ${task.cpus} \\
    $args $partial \\
    --outname ${meta.id} > ${meta.id}_makedb_command_log.txt

    """
}
