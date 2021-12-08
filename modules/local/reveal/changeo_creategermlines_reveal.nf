include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process CHANGEO_CREATEGERMLINES_REVEAL {
    tag "$meta.id"
    label 'process_low'
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
    path(imgt_base) // imgt db

    output:
    tuple val(meta), path("*germ-pass.tsv"), emit: tab
    path("*_command_log.txt"), emit: logs

    script:
    def software = getSoftwareName(task.process)
    """
    CreateGermlines.py -d ${tab} -g dmask \\
    -r ${imgt_base}/${meta.species}/vdj/ --format airr --outdir . \\
    --log ${meta.id}.log --outname ${meta.id} $options.args > "${meta.id}_create-germlines_command_log.txt"
    ParseLog.py -l ${meta.id}.log -f ID V_CALL D_CALL J_CALL
    """
}


