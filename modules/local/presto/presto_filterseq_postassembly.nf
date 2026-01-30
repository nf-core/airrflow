process PRESTO_FILTERSEQ_POSTASSEMBLY {
    tag "$meta.id"
    label "process_medium"
    label 'immcantation'

    conda "bioconda::presto=0.7.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_presto:0f8e73dc0555493d' :
        'biocontainers/presto:0.7.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*quality-pass.fastq") ,  emit: reads
    path "*_command_log.txt" , emit: logs
    path "versions.yml" , emit: versions
    path "*.tab" , emit: log_tab

    script:
    """
    FilterSeq.py quality -s $reads -q ${params.filterseq_q} --outname ${meta.id} --log ${reads.baseName}.log --nproc ${task.cpus} > ${meta.id}_command_log.txt
    ParseLog.py -l ${reads.baseName}.log -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
