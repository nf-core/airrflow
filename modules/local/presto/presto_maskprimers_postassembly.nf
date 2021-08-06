include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_MASKPRIMERS_POSTASSEMBLY {
    tag "$meta.id"
    label "process_high"

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
    path(cprimers)
    path(vprimers)

    output:
    tuple val(meta), path("*REV_primers-pass.fastq") , emit: reads
    path "*command_log.txt", emit: logs
    path "*.log"
    path "*.tab", emit: log_tab


    script:
    def revpr = ''
    if (params.primer_revpr) {
        revpr = '--revpr'
    }
    if (params.cprimer_position == "R1") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${cprimers} --start ${params.cprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > "${meta.id}_command_log.txt"
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${vprimers} --start ${params.vprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> "${meta.id}_command_log.txt"
        ParseLog.py -l "${meta.id}-FWD.log" "${meta.id}-REV.log" -f ID PRIMER ERROR
        """
    } else if (params.cprimer_position == "R2") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${vprimers} --start ${params.cprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > "${meta.id}_command_log.txt"
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${cprimers} --start ${params.vprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> "${meta.id}_command_log.txt"
        ParseLog.py -l "${meta.id}-FWD.log" "${meta.id}-REV.log" -f ID PRIMER ERROR
        """
    } else {
        exit 1, "Error in determining cprimer positon."
    }

}
