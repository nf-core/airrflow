// Import generic module functions
process MERGE_UMI {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.8.0 conda-forge::biopython=1.74"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' :
        'biocontainers/mulled-v2-adc9bb9edc31eb38b3c24786a83b7dfa530e2bea:47d6d7765d7537847ced7dac873190d164146022-0' }"

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
