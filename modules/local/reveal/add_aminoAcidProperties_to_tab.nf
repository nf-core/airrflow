process ADD_AMINOACIDPROPERTIES_TO_TAB {
    tag "$meta.id"
    label 'immcantation'
    label 'process_low'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:3.3.0':
        'docker.io/immcantation/airrflow:3.3.0' }"

    cache 'deep' // Without 'deep' this process would run when using -resume

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*aap-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "versions.yml", emit: versions

    script:
    """
    add_aminoAcidProperties.R --repertoire ${tab} --outname ${meta.id}_aap-pass.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(paste(packageVersion('dplyr'), collapse='.'))")
        optparse: \$(Rscript -e "library(optparse); cat(paste(packageVersion('optparse'), collapse='.'))")
        alakazam: \$(Rscript -e "library(alakazam); cat(paste(packageVersion('alakazam'), collapse='.'))")
        stringr: \$(Rscript -e "library(stringr); cat(paste(packageVersion('stringr'), collapse='.'))")
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
    END_VERSIONS
    """
}
