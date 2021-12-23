process PRESTO_PARSEHEADERS_METADATA {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads

    script:
    """
    ParseHeaders.py add -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f sample_id subject_id species pcr_target_locus -u ${meta.id} ${meta.subject} ${meta.species} ${meta.locus}
    """
}
