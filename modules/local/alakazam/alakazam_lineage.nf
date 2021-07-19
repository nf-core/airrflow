include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process ALAKAZAM_LINEAGE {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-alakazam=1.0.2 bioconda::changeo=1.0.2 bioconda::phylip=3.697" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:ec93fe5ff5457014204d1537f8b85458056510bb-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:ec93fe5ff5457014204d1537f8b85458056510bb-0"                        // Docker image
    }

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format

    output:
    tuple val(meta), path("${tab}"), emit: tab
    path("*.version.txt"), emit: version
    path("*.tsv")
    path("Clone_tree_plots/*.pdf")
    path("Graphml_trees/All_graphs_patient.graphml")

    script:
    def software = getSoftwareName(task.process)
    """
    which dnapars > dnapars_exec.txt
    lineage_reconstruction.R ${tab} $options.args
    merge_graphs.sh
    Rscript -e "library(alakazam); write(x=as.character(packageVersion('alakazam')), file='${software}.version.txt')"
    echo \$(R --version 2>&1) | awk -F' '  '{print \$3}' > R.version.txt
    """

}
