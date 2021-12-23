process CHANGEO_ASSIGNGENES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta

    output:
    path("*igblast.fmt7"), emit: blast
    tuple val(meta), path("$reads"), emit: fasta
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    AssignGenes.py igblast -s $reads -b $igblast --organism $meta.species --loci ${meta.locus.toLowerCase()} $args --nproc $task.cpus --outname "$meta.id"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
