include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_COLLAPSESEQ {
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
    tuple val(meta), path("*_collapse-unique.fastq") , emit: reads
    path("*_command_log.txt") , emit: logs
    path("*.log")
    path("*_table.tab")
    

    script:
    """
    CollapseSeq.py -s $reads -n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname ${meta.id} --log ${meta.id}.log > "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}.log" -f HEADER DUPCOUNT
    """
}