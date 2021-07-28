include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_BUILDCONSENSUS {
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
    tuple val(meta), path("*_R1_consensus-pass.fastq"), path("*_R2_consensus-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*_R1.log")
    path("*_R2.log")
    path("*.tab"), emit: log_tab


    script:
    """
    BuildConsensus.py -s $R1 --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror 0.1 --maxgap 0.5 --outname ${meta.id}_R1 --log ${meta.id}_R1.log > "${meta.id}_command_log.txt"
    BuildConsensus.py -s $R2 --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror 0.1 --maxgap 0.5 --outname ${meta.id}_R2 --log ${meta.id}_R2.log >> "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}_R1.log" "${meta.id}_R2.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
    """
}
