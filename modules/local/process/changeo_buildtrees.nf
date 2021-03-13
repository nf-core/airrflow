include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CHANGEO_BUILDTREES {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge:r-alakazam=1.0.2 bioconda::changeo=1.0.2 bioconda::igphyml=1.1.3" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-d432bd3f78aaba1be2f7eb105c18998acb64d739:65ec671cb141d3fa9117e0c37ac9dcb83970c883-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-d432bd3f78aaba1be2f7eb105c18998acb64d739:65ec671cb141d3fa9117e0c37ac9dcb83970c883-0" // Docker image
    }

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format

    output:
    tuple val(meta), path("*_lineages.tsv")

    script:
    def software = getSoftwareName(task.process)
    """
    BuildTrees.py -d ${tab} --outname ${meta.id} --log ${meta.id}.log --collapse --nproc $task.cpus
    """
}
