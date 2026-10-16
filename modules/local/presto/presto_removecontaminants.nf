process PRESTO_REMOVECONTAMINANTS {
    tag "all_reads_removecontam"
    label "process_medium"
    label 'immcantation'

    //conda "bioconda::presto=0.7.8"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
    //    'biocontainers/presto:0.7.8--pyhdfd78af_0' }"
    container "docker.io/ggabernet/presto:contaminants"

    input:
    path(reads)

    output:
    path("*R1_pair-pass_contam-pass.fastq"), emit: reads_R1
    path("*R2_pair-pass_contam-pass.fastq"), emit: reads_R2
    path "*_command_log.txt", emit: logs
    tuple val("${task.process}"), val('presto'), eval('RemoveContaminants.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args?: ''
    """
    RemoveContaminants.py \
    -1 *R1_pair-pass.fastq \
    -2 *R2_pair-pass.fastq \
    --nproc ${task.cpus} \
    ${args} > remove_contaminants_command_log.txt
    """
}
