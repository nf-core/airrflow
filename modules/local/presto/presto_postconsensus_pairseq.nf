process PRESTO_POSTCONSENSUS_PAIRSEQ {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path("${meta.id}_R1.fastq"), path("${meta.id}_R2.fastq")

    output:
    tuple val(meta), path("*R1_pair-pass.fastq"), path("*R2_pair-pass.fastq") , emit: reads
    path "*_command_log.txt", emit: logs
    tuple val("${task.process}"), val('presto'), eval('PairSeq.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    """
    PairSeq.py -1 ${meta.id}_R1.fastq -2 ${meta.id}_R2.fastq --coord presto > ${meta.id}_command_log.txt

    """
}
