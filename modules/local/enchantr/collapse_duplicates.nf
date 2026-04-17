process COLLAPSE_DUPLICATES {
    tag "$meta.id"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "docker.io/immcantation/airrflow:5.0.0dev"

    input:
    tuple val(meta), path(tabs) // tuple [val(meta), sequence tsv in AIRR format ]

    output:
    tuple val(meta), path("*/*/*collapse-pass.tsv"), emit: tab // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml" , emit: versions

    script:
    """
    Rscript ${projectDir}/bin/reveal_collapseDuplicates.R \\
        --repertoire ${tabs.join(',')} \\
        --collapseby ${params.collapseby} \\
        --ids ${meta.id} \\
        --outname ${meta.id} \\
        > ${meta.id}_collapse_command_log.txt

    mkdir -p ${meta.id}_collapse_report/repertoires
    mv *collapse-pass.tsv ${meta.id}_collapse_report/repertoires/

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
