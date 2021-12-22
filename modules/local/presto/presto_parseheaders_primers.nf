process PRESTO_PARSEHEADERS_PRIMERS {
    tag "$meta.id"
    label "process_low"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::presto=0.6.2=py_0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.6.2--py_0' :
        'quay.io/biocontainers/presto:0.6.2--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads

    script:
    if (params.cprimer_position == "R1") {
        """
        ParseHeaders.py copy -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f $options.args --act first last -k C_PRIMER V_PRIMER
        """
    } else if (params.cprimer_position == "R2") {
        """
        ParseHeaders.py copy -s $reads -o "${reads.baseName}_reheader-pass.fastq" -f $options.args --act first last -k V_PRIMER C_PRIMER
        """
    }

}
