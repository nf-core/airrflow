process REMOVE_CHIMERIC {
    tag "$meta.id"

    label 'process_high'
    label 'process_long_parallelized'

    label 'enchantr'


    conda (params.enable_conda ? "bioconda::r-enchantr=0.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.0.1--r41hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.0.1--r41hdfd78af_0' }"


    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(imgt_base)

    output:
    tuple val(meta), path("*chimera-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" //, emit: chimera_report
    path "versions.yml" , emit: versions

    script:
    """
    Rscript -e "enchantr:::enchantr_report('chimera_analysis', \\
        report_params=list('input'='${tab}',\\
        'outdir'=getwd(), \\
        'nproc'=${task.cpus},\\
        'outname'='${meta.id}', \\
        'log'='${meta.id}_chimeric_command_log'))"

    echo "\"${task.process}\":" > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml

    mv enchantr ${meta.id}_chimera_report
    """
}
