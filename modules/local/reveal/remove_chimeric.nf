include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process REMOVE_CHIMERIC {
    tag "$meta.id"
    label 'immcantation'

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

    script:
    germline_db = tab.getBaseName().toString() + '_germ-pass.tsv'
    """
    CreateGermlines.py -d $tab -r ${imgt_base}/${meta.species}/vdj/ -g dmask --format airr > "${meta.id}_${task.process}_create-germlines_command_log.txt"
    reveal_chimeric.R --repertoire ${germline_db} --outname ${meta.id} > "${meta.id}_${task.process}_chimeric_command_log.txt"
    """
}
