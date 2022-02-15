process CHANGEO_ASSIGNGENES_REVEAL {
    tag "$meta.id"
    label 'process_high'
    label 'immcantation'

    // TODO: update container
    container "immcantation/suite:devel"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta

    output:
    path("*igblast.fmt7"), emit: blast
    tuple val(meta), path("$reads"), emit: fasta
    path "versions.yml" , emit: versions
    path "*_command_log.txt" , emit: logs

    script:
    def args = task.ext.args ?: ''
    """
    AssignGenes.py igblast -s $reads -b $igblast --organism "$meta.species" --loci "$meta.locus" --format blast --nproc $task.cpus --outname "$meta.id"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS"""
}
