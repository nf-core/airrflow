include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_PARSEHEADERS_PRIMERS {
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
    def field_names = "PRIMER PRIMER"
    if (params.umi) {
        field_names = "PRCONS PRCONS"
    }
    if (params.cprimer_position == "R1") {
        """
        ParseHeaders.py copy -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f $field_names --act first last -k C_PRIMER V_PRIMER
        """
    } else if (params.cprimer_position == "R2") {
        """
        ParseHeaders.py copy -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f $field_names --act first last -k V_PRIMER C_PRIMER
        """
    }

}
