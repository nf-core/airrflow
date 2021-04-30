include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_PAIRSEQ {
    tag "$meta.id"
    label "process_low"

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
    tuple val(meta), path("${meta.id}_R1.fastq"), path("${meta.id}_R2.fastq")

    output:
    tuple val(meta), path("*R1_pair-pass.fastq"), path("*R2_pair-pass.fastq") , emit: reads
    path "*_command_log.txt", emit: logs

    script:
    def copyfield = (params.index_file & params.umi_position == 'R1') ? "--1f BARCODE" : "--2f BARCODE"
    """
    PairSeq.py -1 '${meta.id}_R1.fastq' -2 '${meta.id}_R2.fastq' $copyfield --coord illumina > "${meta.id}_command_log.txt"
    """
}