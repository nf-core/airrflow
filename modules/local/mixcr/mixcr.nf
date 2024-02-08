process MIXCR_MIXCR {
    tag "$meta.id"
    label 'process_medium'

    secret 'MIXCR_LICENSE'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO check if you can use their conda for conda
    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0' }"

    input:
    tuple val(meta), path(reads)
    path(imgt_json) // imgt db
    val(kit)

    output:
    tuple val(meta), path('*.json'), emit: json
    tuple val(meta), path('*.clns'), emit: clns
    tuple val(meta), path('*.vdjca'), emit: vdjca
    tuple val(meta), path('*.txt'), emit: txt
    tuple val(meta), path('*.tsv'), emit: tsv
    tuple val(meta), path('*clones*.tsv'), emit: clones_tsv
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
    """
    # activate license
    if [ \${MIXCR_LICENSE:-"unset"} != "unset" ]; then
        echo "Initializing MIXCR_LICENSE env variable"
        export MI_LICENSE=\$MIXCR_LICENSE
    fi

    mixcr analyze ${kit} \\
        --library imgt \\
        ${reads} \\
        ${prefix} \\
        -t ${task.cpus} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version |& sed '1!d ; s/mixcr //')
    END_VERSIONS
    """
}
