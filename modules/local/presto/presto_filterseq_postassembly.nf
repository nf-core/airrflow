include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)


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
    FilterSeq.py quality -s $reads -q ${params.filterseq_q} --outname "${meta.id}" --log "${reads.baseName}.log" --nproc ${task.cpus} > "${meta.id}_command_log.txt"
    ParseLog.py -l "${reads.baseName}.log" -f ID QUALITY
    FilterSeq.py --version | awk -F' '  '{print \$2}' > ${software}.version.txt
    """
}
