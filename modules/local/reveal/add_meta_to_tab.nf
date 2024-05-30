process ADD_META_TO_TAB {
    tag "$meta.id"
    label 'immcantation'
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/immcantation/airrflow:4.1.0':
        'docker.io/immcantation/airrflow:4.1.0' }"

    cache 'deep' // Without 'deep' this process would run when using -resume

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format
    path(validated_input)

    output:
    tuple val(meta), path("*meta-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "versions.yml", emit: versions

    script:
    """
    reveal_add_metadata.R --repertoire ${tab} --metadata ${validated_input} --input_id ${meta.id} --outname ${meta.id} > ${meta.id}_add-meta_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(paste(packageVersion('dplyr'), collapse='.'))")
        optparse: \$(Rscript -e "library(optparse); cat(paste(packageVersion('optparse'), collapse='.'))")
        R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
    END_VERSIONS
    """
}
