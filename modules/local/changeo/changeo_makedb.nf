process CHANGEO_MAKEDB {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'


    conda "bioconda::changeo=1.3.5 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0 bioconda::airr=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8d583137003a25ebc5fdd51c08fccacc3f3673c9fc162c92b006335a2cd6df0/data' :
        'community.wave.seqera.io/library/airr_changeo_igblast_wget:7303fb2eff967b2e' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect()
    path(reference_fasta)

    output:
    tuple val(meta), path("*db-pass.tsv.gz"), emit: tab //compressed sequence table in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('changeo'), eval('MakeDb.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions

    script:
    def args = task.ext.args ?: ''
    def partial = meta.species.toLowerCase()=='mouse' && meta.locus.toLowerCase()=='tr'  ? '--partial' : ''
    """
    MakeDb.py igblast -i $igblast -s $reads -r \\
    ${reference_fasta}/${meta.species.toLowerCase()}/vdj/ \\
    --nproc ${task.cpus} --gzip-output \\
    $args $partial \\
    --outname ${meta.id} > ${meta.id}_makedb_command_log.txt

    """
}
