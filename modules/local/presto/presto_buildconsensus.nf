process PRESTO_BUILDCONSENSUS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)
    val buildconsensus_maxerror
    val buildconsensus_maxgap
    val use_consensus_R1
    val use_consensus_R2
    val primer_consensus
    val cluster_sets

    output:
    tuple val(meta), path("*_R1_consensus-pass.fastq"), path("*_R2_consensus-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*.tab"), emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('BuildConsensus.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def cluster_sets_val = cluster_sets ? "--bf CLUSTER" : "--bf BARCODE"
    def primer_consensus_R1 = use_consensus_R1 ? "--prcons $primer_consensus" : ""
    def primer_consensus_R2 = use_consensus_R2 ? "--prcons $primer_consensus" : ""
    """
    BuildConsensus.py -s $R1 \\
    --nproc ${task.cpus} \\
    --maxerror ${buildconsensus_maxerror} \\
    --maxgap ${buildconsensus_maxgap} \\
    ${primer_consensus_R1} \\
    ${cluster_sets_val} \\
    ${args} \\
    --outname ${meta.id}_R1 \\
    --log ${meta.id}_R1.log > ${meta.id}_command_log.txt
    BuildConsensus.py -s $R2 --nproc ${task.cpus} \\
    --maxerror ${buildconsensus_maxerror} \\
    --maxgap ${buildconsensus_maxgap} \\
    ${primer_consensus_R2} \\
    ${cluster_sets_val} \\
    ${args2} \\
    --outname ${meta.id}_R2 \\
    --log ${meta.id}_R2.log >> ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log ${args3}

    """
}
