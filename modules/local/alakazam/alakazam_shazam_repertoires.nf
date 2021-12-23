process ALAKAZAM_SHAZAM_REPERTOIRES {
    tag "report"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-alakazam=1.0.2 conda-forge::r-shazam=1.0.2 conda-forge::r-kableextra=1.3.4 conda-forge::r-knitr=1.33 conda-forge::r-stringr=1.4.0 conda-forge::r-dplyr=1.0.6" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3420a264d7f8006cc73fc3c3843d4545b235404f:4cc818718337222966d00ce39968abea1c328367-0' :
        'quay.io/biocontainers/mulled-v2-3420a264d7f8006cc73fc3c3843d4545b235404f:4cc818718337222966d00ce39968abea1c328367-0' }"

    input:
    path(tab) // sequence tsv table in AIRR format
    path("Table_sequences.tsv")
    path(repertoire_report)

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
