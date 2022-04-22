process PRESTO_PAIRSEQ {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::presto=0.7.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.0--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path("${meta.id}_R1.fastq"), path("${meta.id}_R2.fastq")

    output:
    tuple val(meta), path("*R1_pair-pass.fastq"), path("*R2_pair-pass.fastq") , emit: reads
    path "*_command_log.txt", emit: logs
    path "versions.yml" , emit: versions

    script:
    def copyfield = (params.index_file | params.umi_position == 'R1') ? "--1f BARCODE" : "--2f BARCODE"
    """
    PairSeq.py -1 '${meta.id}_R1.fastq' -2 '${meta.id}_R2.fastq' $copyfield --coord illumina > "${meta.id}_command_log.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( PairSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
