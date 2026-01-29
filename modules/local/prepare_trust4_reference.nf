process PREPARE_TRUST4_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(R1), path(R2)
    path(reference_igblast)

    output:
    path("trust4_reference.fa") , emit: trust4_reference
    path "versions.yml"         , emit: versions

    script:
    """
    cat ${reference_igblast}/fasta/imgt_${meta.species.toLowerCase()}_*.fasta \\
    ${reference_igblast}/fasta/imgt_${meta.species.toLowerCase()}_*.fasta >> trust4_reference.fa

    cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
    """


}
