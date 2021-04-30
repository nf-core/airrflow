include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_PARSEHEADERS_METADATA {
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
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads
    
    script:
    """
    ParseHeaders.py add -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f SAMPLE_CODE SOURCE TREATMENT EXTRACT_TIME POPULATION -u ${meta.id} ${meta.source} ${meta.treatment} ${meta.time} ${meta.population}
    """
}