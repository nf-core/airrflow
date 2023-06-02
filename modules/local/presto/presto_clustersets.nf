process PRESTO_CLUSTERSETS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*_R1_cluster-pass.fastq"), path("*_R2_cluster-pass.fastq"), emit: reads
    path "*_command_log.txt", emit: logs
    path "*.log"
    path "*.tab", emit: log_tab
    path("versions.yml"), emit: versions

    script:
    """
    ClusterSets.py set --nproc ${task.cpus} -s $R1 --outname ${meta.id}_R1 --exec vsearch --log ${meta.id}_R1.log > ${meta.id}_command_log.txt
    ClusterSets.py set --nproc ${task.cpus} -s $R2 --outname ${meta.id}_R2 --exec vsearch --log ${meta.id}_R2.log >> ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log -f ID BARCODE SEQCOUNT CLUSTERS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( ClusterSets.py --version | awk -F' '  '{print \$2}' )
        vsearch: \$( vsearch --version &> vsearch.txt; cat vsearch.txt | head -n 1 | grep -o 'v[0-9\\.]\\+' )
    END_VERSIONS
    """
}
