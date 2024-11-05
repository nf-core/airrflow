process CHANGEO_ASSIGNGENES {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.0 bioconda::igblast=1.22.0 conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' :
        'biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast references

    output:
    path("*igblast.fmt7"), emit: blast
    tuple val(meta), path("$reads"), emit: fasta
    path "versions.yml" , emit: versions
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def args = task.ext.args ?: ''
    """
    AssignGenes.py igblast \\
    -s $reads \\
    -b ${igblast} \\
    --vdb "${meta.reference}_${meta.species}_${meta.locus.toLowerCase()}_v" \\
    --ddb "${meta.reference}_${meta.species}_${meta.locus.toLowerCase()}_d" \\
    --jdb "${meta.reference}_${meta.species}_${meta.locus.toLowerCase()}_j" \\
    --cdb "${meta.reference}_${meta.species}_${meta.locus.toLowerCase()}_c" \\
    --organism $meta.species \\
    --loci ${meta.locus.toLowerCase()} \\
    $args \\
    --nproc $task.cpus \\
    --outname $meta.id > ${meta.id}_changeo_assigngenes_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
