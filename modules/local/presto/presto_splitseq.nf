include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_SPLITSEQ {
    tag "$meta.id"

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
    tuple val(meta), path("*_atleast-2.fasta"), emit: fasta
    path("*_command_log.txt"), emit: logs

    script:
    def field_option = "-f DUPCOUNT"
    if (params.umi) {
        field_option = "-f CONSCOUNT"
    }
    """
    SplitSeq.py group -s $reads --outname ${meta.id} $options.args --fasta > "${meta.id}_command_log.txt"
    """
}
