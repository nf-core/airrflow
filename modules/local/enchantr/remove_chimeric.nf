include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process REMOVE_CHIMERIC {
    tag "$meta.id"
    label 'immcantation'
    label 'enchantr'
    label 'process_high'
    label 'process_long'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
    }

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(imgt_base)

    output:
    tuple val(meta), path("*chimera-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report", emit: chimera_report

    script:
    germline_db = tab.getBaseName().toString() + '_germ-pass.tsv'
    """
    CreateGermlines.py -d "${tab}" -r ${imgt_base}/${meta.species}/vdj/ -g dmask --format airr --outdir . > "${meta.id}_create-germlines_command_log.txt"
    Rscript -e "enchantr:::enchantr_report('chimera_analysis', report_params=list('input'='${germline_db}','outdir'=getwd(), 'nproc'=${task.cpus},'outname'='${meta.id}', 'log'='${meta.id}_chimeric_command_log'))"
    mv enchantr ${meta.id}_chimera_report
    """
}
