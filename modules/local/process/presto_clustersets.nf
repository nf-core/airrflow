include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_CLUSTERSETS {
    tag "$meta.id"
    label "process_long_parallelized"

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
    tuple val(meta), path("*_R1_cluster-pass.fastq"), path("*_R2_cluster-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    //path("*.version.txt"), emit: version

    script:
    """
    ClusterSets.py set --nproc ${task.cpus} -s $R1 --outname ${meta.id}_R1 --exec vsearch > "${meta.id}_command_log.txt"
    ClusterSets.py set --nproc ${task.cpus} -s $R2 --outname ${meta.id}_R2 --exec vsearch >> "${meta.id}_command_log.txt"
    """
    //TODO add version scraping for vsearch
    // tried with:  vsearch --version > vsearch.txt; cat vsearch.txt | head -n 1 | grep -Eo "v[0-9\.]{4,7}" > vsearch.version.txt

}