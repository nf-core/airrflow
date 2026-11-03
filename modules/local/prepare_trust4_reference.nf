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
    tuple val("${task.process}"), val('cat'), eval('cat --version 2>&1 | head -n 1 | sed \'s/^.*coreutils) //; s/ .*$//\''), emit: versions_cat, topic: versions

    script:
    """
    cat ${reference_igblast}/fasta/${meta.species.toLowerCase()}_*.fasta >> trust4_reference.fa

    """


}
