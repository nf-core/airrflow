process MIXCR_OVERLAP_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    secret 'MIXCR_LICENSE'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0' }"

    input:
    tuple val(meta), path(mixcr_overlap_json)

    output:
    tuple val(meta), path('*')  , emit: outs

    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MiXCR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # activate license
    if [ \${MIXCR_LICENSE:-"unset"} != "unset" ]; then
        echo "Initializing MIXCR_LICENSE env variable"
        export MI_LICENSE=\$MIXCR_LICENSE
    fi

    mixcr exportPlots overlap \\
        ${mixcr_overlap_json} \\
        ${prefix}.overlap.pdf




    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version |& sed '1!d ; s/mixcr //')
    END_VERSIONS
    """
}
