process ALAKAZAM_SHAZAM_REPERTOIRES {
    tag "report"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-base=4.1.2 bioconda::r-alakazam=1.2.0 bioconda::r-shazam=1.1.0 conda-forge::r-kableextra=1.3.4 conda-forge::r-knitr=1.33 conda-forge::r-stringr=1.4.0 conda-forge::r-dplyr=1.0.6 conda-forge::r-optparse=1.7.1" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7da73314bcc47157b442d16c3dcfbe81e75a404f:9bb35f8114dffcd97b3afb5de8587355aca16b66-0' :
        'quay.io/biocontainers/mulled-v2-7da73314bcc47157b442d16c3dcfbe81e75a404f:9bb35f8114dffcd97b3afb5de8587355aca16b66-0' }"

    input:
    path(tab) // sequence tsv table in AIRR format
    path("Table_sequences.tsv")
    tuple path(repertoire_report), path(references), path(css), path(logo)

    output:
    path("versions.yml"), emit: versions
    path("repertoire_comparison")
    path("Bcellmagic_report.html")

    script:
    """
    execute_report.R ${repertoire_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        shazam: \$(Rscript -e "library(shazam); cat(paste(packageVersion('shazam'), collapse='.'))")
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
    END_VERSIONS
    """
}
