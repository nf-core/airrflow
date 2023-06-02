process ADD_META_TO_TAB {
    tag "$meta.id"
    label 'immcantation'
    label 'process_single'

    conda "bioconda::r-enchantr=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-enchantr:0.1.1--r42hdfd78af_0':
        'quay.io/biocontainers/r-enchantr:0.1.1--r42hdfd78af_0' }"

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
