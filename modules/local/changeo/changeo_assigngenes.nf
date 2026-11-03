process CHANGEO_ASSIGNGENES {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast references

    output:
    path("*igblast.fmt7"), emit: blast
    tuple val(meta), path("$reads"), emit: fasta
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('changeo'), eval('AssignGenes.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def args = task.ext.args ?: ''
    def species = meta.species.toLowerCase()
    def loci = meta.locus.toLowerCase()
    """
    AssignGenes.py igblast \\
        -s $reads \\
        -b $igblast \\
        --organism $meta.species \\
        --loci ${loci} \\
        --vdb ${species}_${loci}_v \\
        --ddb ${species}_${loci}_d \\
        --jdb ${species}_${loci}_j \\
        --cdb ${species}_${loci}_c \\
        $args \\
        --nproc $task.cpus \\
        --outname $meta.id > ${meta.id}_changeo_assigngenes_command_log.txt

    """
}
