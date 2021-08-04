include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label "process_medium"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0"  // Singularity image
    } else {
        container "quay.io/biocontainers/presto:0.6.2--py_0"                        // Docker image
    }

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*R1_quality-pass.fastq"), path("*R2_quality-pass.fastq") ,  emit: reads
    path "*_command_log.txt" , emit: logs
    path "*.version.txt" , emit: version
    path "*_R1.log"
    path "*_R2.log"
    path "*.tab" , emit: log_tab

    script:
    def software = getSoftwareName(task.process)
    """
    FilterSeq.py quality -s $R1 -q ${params.filterseq_q} --outname "${meta.id}_R1" --log "${R1.baseName}_R1.log" --nproc ${task.cpus} > "${meta.id}_command_log.txt"
    FilterSeq.py quality -s $R2 -q ${params.filterseq_q} --outname "${meta.id}_R2" --log "${R2.baseName}_R2.log" --nproc ${task.cpus} >> "${meta.id}_command_log.txt"
    ParseLog.py -l "${R1.baseName}_R1.log" "${R2.baseName}_R2.log" -f ID QUALITY
    FilterSeq.py --version | awk -F' '  '{print \$2}' > ${software}.version.txt
    """
}

// Filter single end reads if after assembly (will be the typical case without UMIs)
process PRESTO_FILTERSEQ_POSTASSEMBLY {
    tag "$meta.id"
    label "process_medium"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0"  // Singularity image
    } else {
        container "quay.io/biocontainers/presto:0.6.2--py_0"                        // Docker image
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*quality-pass.fastq") ,  emit: reads
    path "*_command_log.txt" , emit: logs
    path "*.version.txt" , emit: version
    path "*.log"
    path "*.tab" , emit: log_tab

    script:
    def software = getSoftwareName(task.process)
    """
    FilterSeq.py quality -s $R1 -q ${params.filterseq_q} --outname "${meta.id}" --log "${R1.baseName}.log" --nproc ${task.cpus} > "${meta.id}_command_log.txt"
    ParseLog.py -l "${R1.baseName}.log" -f ID QUALITY
    FilterSeq.py --version | awk -F' '  '{print \$2}' > ${software}.version.txt
    """
}
