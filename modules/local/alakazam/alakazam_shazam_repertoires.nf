include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process ALAKAZAM_SHAZAM_REPERTOIRES {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-alakazam=1.0.2 conda-forge::r-shazam=1.0.2" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/quay.io/biocontainers/mulled-v2-62e1fba87191802d6c6cd9a0faee2e3aba60bcae:bd564bd32d225d1eab3f128ae5276187ea6730d1-0"  // Singularity image
    } else {
        container "nfcore/bcellmagic:2.0dev"                        // Docker image
    }

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    path("Table_sequences.tsv")

    output:
    tuple val(meta), path("${tab}"), emit: tab
    path("*.version.txt"), emit: version
    path("repertoire_comparison/*")
    path("Bcellmagic_report.html")

    script:
    def software = getSoftwareName(task.process)
    """
    execute_report.R "$projectDir/assets/repertoire_comparison.Rmd"
    Rscript -e "library(alakazam); write(x=as.character(packageVersion('alakazam')), file='${software}.version.txt')"
    Rscript -e "library(shazam); write(x=as.character(packageVersion('shazam')), file='shazam.version.txt')"
    echo \$(R --version 2>&1) | awk -F' '  '{print \$3}' > R.version.txt
    """
}