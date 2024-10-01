process MIXCR_OVERLAP_POSTANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    secret 'MIXCR_LICENSE'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0' }"

    input:
    tuple val(meta), path(clns)
    val(downsampling)
    val(weight_function)
    val(productive)
    val(drop_outliers)
    val(criteria)
    path(imgt_json) // imgt db


    output:
    tuple val(meta), path('*')      , emit: outs
    tuple val(meta), path('*.json') , emit: mixcr_overlap_json

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MiXCR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def productive_only = productive ? '--only-productive' : ''
    def drop_outliers = drop_outliers ? '--drop-outliers' : ''
    def criteria_pick = criteria ? "--criteria '${criteria}'" : ''
    """
    # activate license
    if [ \${MIXCR_LICENSE:-"unset"} != "unset" ]; then
        echo "Initializing MIXCR_LICENSE env variable"
        export MI_LICENSE=\$MIXCR_LICENSE
    fi

    mixcr postanalysis overlap \\
        --default-downsampling ${downsampling} \\
        --default-weight-function ${weight_function} \\
        ${productive_only} \\
        ${drop_outliers} \\
        ${criteria_pick} \\
        ${clns} \\
        ${prefix}.overlap_postanalysis.json \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version |& sed '1!d ; s/mixcr //')
    END_VERSIONS
    """
}
