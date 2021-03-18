include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SHAZAM_TIGGER_THRESHOLD {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-shazam=1.0.2 conda-forge::r-tigger=1.0.0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-02abbf3250c7d008f9d739d7e72ff16a64ae95fc:a440898f307c775f85c44c2812cb28afd74f011d-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-02abbf3250c7d008f9d739d7e72ff16a64ae95fc:a440898f307c775f85c44c2812cb28afd74f011d-0"                        // Docker image
    }

    input:
    tuple val(meta), path(tab) // tsv tab in AIRR format
    path(imgt_base) // igblast fasta

    output:
    tuple val(meta), path("*genotyped.tab"), emit: tab
    path("threshold.txt"), emit: threshold
    path("*.version.txt") , emit: version
    path("*genotype.fasta"), emit: fasta
    path("*/*.pdf")

    script:
    def software = getSoftwareName(task.process)
    if (params.loci == 'ig'){
        """
        TIgGER-shazam.R $tab $params.loci $params.threshold_method ${imgt_base}/${params.species}/vdj/imgt_human_IGHV.fasta 
        """
    } else if (params.loci == 'tr'){
        """
        TIgGER-shazam.R $tab $params.loci $params.threshold_method \\
        "${imgt_base}/${params.species}/vdj/imgt_${params.species}_TRAV.fasta" \\
        "${imgt_base}/${params.species}/vdj/imgt_${params.species}_TRBV.fasta" \\
        "${imgt_base}/${params.species}/vdj/imgt_${params.species}_TRDV.fasta" \\
        "${imgt_base}/${params.species}/vdj/imgt_${params.species}_TRGV.fasta"
        Rscript -e "library(shazam); write(x=as.character(packageVersion('shazam')), file='${software}.version.txt')"
        Rscript -e "library(tigger); write(x=as.character(packageVersion('tigger')), file='tigger.version.txt')"
        """
    }
}
