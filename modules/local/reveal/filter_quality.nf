process FILTER_QUALITY {
    tag "$meta.id"
    label 'immcantation'
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:3.2.0':
        'docker.io/immcantation/airrflow:3.2.0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*quality-pass.tsv"), optional:true, emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "versions.yml", emit: versions

    script:
    """
    reveal_filter_quality.R --repertoire $tab --outname ${meta.id} > ${meta.id}_fq_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        optparse: \$(Rscript -e "library(optparse); cat(paste(packageVersion('optparse'), collapse='.'))")
        stringi: \$(Rscript -e "library(stringi); cat(paste(packageVersion('stringi'), collapse='.'))")
        dplyr: \$(Rscript -e "library(dplyr); cat(paste(packageVersion('dplyr'), collapse='.'))")
        airr: \$(Rscript -e "library(airr); cat(paste(packageVersion('airr'), collapse='.'))")
        DT: \$(Rscript -e "library(DT); cat(paste(packageVersion('DT'), collapse='.'))")
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
    END_VERSIONS
    """
}
