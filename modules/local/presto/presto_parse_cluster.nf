include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_PARSE_CLUSTER {
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
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*R1_cluster-pass_reheader.fastq"), path("*R2_cluster-pass_reheader.fastq"), emit: reads
    path("*_log.txt"), emit: logs

    script:
    """
    ParseHeaders.py copy -s $R1 -f BARCODE -k CLUSTER --act cat > "${meta.id}_command_log.txt"
    ParseHeaders.py copy -s $R2 -f BARCODE -k CLUSTER --act cat >> "${meta.id}_command_log.txt"
    """
}