process AIRRFLOW_REPORT {
    tag "${meta.id}"
    label 'process_high'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:4.4.0':
        'docker.io/immcantation/airrflow:4.4.0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    path("Table_sequences.tsv")
    path("Table_sequences_assembled.tsv")
    path(repertoire_report)
    path(css)
    path(logo)

    output:
    path "versions.yml" , emit: versions
    path("repertoire_comparison"), emit: results_folder
    path("*.html"), emit: report_html

    script:
    """
    execute_report.R --report_file ${repertoire_report}

    mkdir repertoire_comparison/repertoires
    cp *clone-pass.tsv repertoire_comparison/repertoires/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        shazam: \$(Rscript -e "library(shazam); cat(paste(packageVersion('shazam'), collapse='.'))")
        stringr: \$(Rscript -e "library(stringr); cat(paste(packageVersion('stringr'), collapse='.'))")
        dplyr: \$(Rscript -e "library(dplyr); cat(paste(packageVersion('dplyr'), collapse='.'))")
        knitr: \$(Rscript -e "library(knitr); cat(paste(packageVersion('knitr'), collapse='.'))")
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
    END_VERSIONS
    """
}
