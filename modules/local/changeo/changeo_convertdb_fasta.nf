include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process CHANGEO_CONVERTDB_FASTA {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    if (!params.custom_container) {
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
        } else {
            container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
        }
    } else {
        container params.custom_container
    }
    
    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*.fasta"), emit: fasta // sequence tsv in AIRR format
    path "*.version.txt" , emit: version
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def software = getSoftwareName(task.process)
    """
    ConvertDb.py fasta -d $tab $options.args --outname ${meta.id}  > "${meta.id}_${task.process}_command_log.txt"
    ConvertDb.py --version | awk -F' '  '{print \$2}' > ${software}.version.txt
    """
}
