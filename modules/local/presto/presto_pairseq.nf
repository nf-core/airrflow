process PRESTO_PAIRSEQ {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_presto:0f8e73dc0555493d' :
        'biocontainers/presto:0.7.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path("${meta.id}_R1.fastq"), path("${meta.id}_R2.fastq")
    val(barcode_position)

    output:
    tuple val(meta), path("*R1_pair-pass.fastq"), path("*R2_pair-pass.fastq") , emit: reads
    path "*_command_log.txt", emit: logs
    path "versions.yml" , emit: versions

    script:
    def copyfield = (barcode_position == 'R1')? '--1f BARCODE' : (barcode_position == 'R2')? '--2f BARCODE' : (barcode_position == 'R1R2')? '--1f BARCODE --2f BARCODE' : (barcode_position == 'clustersets')? '--1f CLUSTER --2f CLUSTER' : ''
    def args = task.ext.args?: ''
    """
    PairSeq.py -1 ${meta.id}_R1.fastq -2 ${meta.id}_R2.fastq $copyfield $args > ${meta.id}_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( PairSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
