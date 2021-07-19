include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process PRESTO_ASSEMBLEPAIRS {
    tag "$meta.id"
    label 'process_long_parallelized'

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
    tuple val(meta), path("*_assemble-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*.log")
    path("*_table.tab")

    script:
    script_options = ''
    parse_options = ''
    if (params.umi) {
        script_options = '--1f CONSCOUNT PRCONS --2f CONSCOUNT PRCONS'
        parse_options = 'ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT LENGTH OVERLAP ERROR PVALUE'
    } else {
        parse_options = 'ID SEQCOUNT PRIMER PRCOUNT PRFREQ LENGTH OVERLAP ERROR PVALUE'
    }
    """
    AssemblePairs.py align -1 $R1 -2 $R2 --nproc ${task.cpus} --coord ${params.assemblepairs_coord} \
        --rc tail --maxerror ${params.assemblepairs_maxerror} \
        ${script_options} \
        --outname ${meta.id} --log ${meta.id}.log > ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}.log -f ${parse_options}
    """
}
