process ALAKAZAM_LINEAGE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-base=4.1.2 bioconda::r-alakazam=1.2.0 bioconda::changeo=1.2.0 bioconda::phylip=3.697" : null)  // Please also update the phylip version manually in the script section below as phylip does not print the version
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:d6b54b0fae81ec8e55d41b5ea9cc8f39d75cf2d7-0' :
        'quay.io/biocontainers/mulled-v2-afe1e5f3879e265b14ec08dd3a1875df9c23630d:d6b54b0fae81ec8e55d41b5ea9cc8f39d75cf2d7-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format

    output:
    tuple val(meta), path("${tab}"), emit: tab
    path "versions.yml" , emit: versions
    path("*.tsv")
    path("Clone_tree_plots/*.pdf"), emit: graph_plots optional true
    path("Graphml_trees/All_graphs_patient.graphml"), emit: graph_export optional true

    script:
    def args = task.ext.args ?: ''
    """
    which dnapars > dnapars_exec.txt
    lineage_reconstruction.R --repertoire ${tab} $args
    merge_graphs.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
        pyhlip: 3.697
    END_VERSIONS
    """

}
