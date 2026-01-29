process PRESTO_PARSE_CLUSTER {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c538a0e310c303233164cbe486b5e5f6bddcf18975a9b20ac2f590f151f03e62/data' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*R1_cluster-pass_reheader.fastq"), path("*R2_cluster-pass_reheader.fastq"), emit: reads
    path("*_log.txt"), emit: logs
    path "versions.yml" , emit: versions


    script:
    """
    ParseHeaders.py copy -s $R1 -f BARCODE -k CLUSTER --act cat > ${meta.id}_command_log.txt
    ParseHeaders.py copy -s $R2 -f BARCODE -k CLUSTER --act cat >> ${meta.id}_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( ParseHeaders.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
