process PREPARE_TRUST4_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::trust4=1.0.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trust4:1.0.13--h43eeafb_0':
        'biocontainers/trust4:1.0.13--h43eeafb_0' }"

    input:
    tuple val(meta), path(R1), path(R2)
    path(reference_igblast)

    output:
    tuple val(meta), path("trust4_reference.fa") , emit: trust4_reference

    script:
    """
    cat ${reference_igblast}/fasta/imgt_${meta.species.toLowerCase()}_*.fasta \\
    ${reference_igblast}/fasta/imgt_${meta.species.toLowerCase()}_*.fasta >> trust4_reference.fa
    """


}
