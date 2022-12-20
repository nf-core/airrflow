process CHANGEO_BUILDTREES {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'


    conda (params.enable_conda ? "conda-forge::r-base=4.1.2 bioconda:r-alakazam=1.2.0 bioconda::changeo=1.2.0 bioconda::igphyml=1.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d432bd3f78aaba1be2f7eb105c18998acb64d739:2c83ca89e577c8839f746f0fe4a6c63ef5984b99-0' :
        'quay.io/biocontainers/mulled-v2-d432bd3f78aaba1be2f7eb105c18998acb64d739:2c83ca89e577c8839f746f0fe4a6c63ef5984b99-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format

    output:
    tuple val(meta), path("*_lineages.tsv")
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    BuildTrees.py -d ${tab} --outname ${meta.id} --log ${meta.id}.log --nproc $task.cpus $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        changeo: \$(AssignGenes.py --version | awk -F' '  '{print \$2}')
        igphyml: \$(igphyml --version | grep -o "IgPhyML [0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+")
    END_VERSIONS
    """
}
