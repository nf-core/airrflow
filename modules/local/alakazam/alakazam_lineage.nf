process ALAKAZAM_LINEAGE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-alakazam=1.0.2 bioconda::changeo=1.0.2 bioconda::phylip=3.697" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quay.io/biocontainers/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:8db5a7467e440fdacef5c725c4e18726a4f0085c-0' :
        'quay.io/biocontainers/quay.io/biocontainers/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:8db5a7467e440fdacef5c725c4e18726a4f0085c-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format

    output:
    tuple val(meta), path("${tab}"), emit: tab
    path("versions.yml"), emit: versions
    path("*.tsv")
    path("Clone_tree_plots/*.pdf")
    path("Graphml_trees/All_graphs_patient.graphml")

    script:
    def args = task.ext.args ?: ''
    """
    which dnapars > dnapars_exec.txt
    lineage_reconstruction.R ${tab} $args
    merge_graphs.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
    END_VERSIONS
    """

}
