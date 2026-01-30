process PRESTO_CLUSTERSETS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda "bioconda::presto=0.7.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_presto:0f8e73dc0555493d' :
        'biocontainers/presto:0.7.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*_R1_cluster-pass.fastq"), path("*_R2_cluster-pass.fastq"), emit: reads
    path "*_command_log.txt", emit: logs
    path "*.tab", emit: log_tab
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    ClusterSets.py set --nproc ${task.cpus} -s $R1 --outname ${meta.id}_R1 $args --log ${meta.id}_R1.log > ${meta.id}_command_log.txt
    ClusterSets.py set --nproc ${task.cpus} -s $R2 --outname ${meta.id}_R2 $args --log ${meta.id}_R2.log >> ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( ClusterSets.py --version | awk -F' '  '{print \$2}' )
        vsearch: \$( vsearch --version &> vsearch.txt; cat vsearch.txt | head -n 1 | grep -o 'v[0-9\\.]\\+' )
    END_VERSIONS
    """
}
