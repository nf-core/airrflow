include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process ADD_META_TO_TAB {
    tag "$meta.id"
    label 'immcantation'
    label 'single_cpu'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
    }

    cache 'deep' // Without 'deep' this process would run when using -resume

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(validated_input)

    output:
    tuple val(meta), path("*meta-pass.tsv"), emit: tab // sequence tsv in AIRR format

    script:
    """
    # TODO: remove not relevant fields
    reveal_add_metadata.R --repertoire "${tab}" --metadata "${validated_input}" --input_id "${meta.id}" --outname "${meta.id}" > "${meta.id}_${task.process}_command_log.txt"
    """
}
