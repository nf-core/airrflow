// Import generic module functions
process MERGE_UMI {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(R1), path(R2), path(I1)

    output:
    tuple val(meta), path('*_R1.fastq.gz'), path('*_R2.fastq.gz')   , emit: reads
    path "versions.yml" , emit: versions

    script:
    """
    merge_R1_umi.py -R1 "${R1}" -I1 "${I1}" -o UMI_R1.fastq.gz --umi_start $params.umi_start --umi_length $params.umi_length
    mv "UMI_R1.fastq.gz" "${meta.id}_UMI_R1.fastq.gz"
    mv "${R2}" "${meta.id}_R2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$(python --version | grep -o "[0-9\\. ]\\+") )
        biopython: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)"))
    END_VERSIONS
    """
}
