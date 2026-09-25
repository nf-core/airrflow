process REMOVE_CHIMERIC {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.0dev"


    input:
    tuple val(meta), path(tab) // compressed sequence tsv in AIRR format
    path(reference_fasta)

    output:
    tuple val(meta), path("*chimera-pass.tsv.gz"), emit: tab // compressed sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report" //, emit: chimera_report
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    """
    Rscript -e "enchantr:::enchantr_report('chimera_analysis', \\
        report_params=list('input'='${tab}',\\
        'outdir'=getwd(), \\
        'nproc'=${task.cpus},\\
        'outname'='${meta.id}', \\
        'log'='${meta.id}_chimeric_command_log'))"

    cp -r enchantr ${meta.id}_chimera_report && rm -rf enchantr

    """
}
