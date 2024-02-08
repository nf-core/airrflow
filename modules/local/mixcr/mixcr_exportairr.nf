process MIXCR_MIXCREXPORTAIRR {
    tag "$meta.id"
    label 'process_single'

    secret 'MIXCR_LICENSE'


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0' }"

    input:

    tuple val(meta), path(clns)
    path(imgt_json) // imgt db
    
    output:
    tuple val(meta), path("*.airr.tsv"), emit: mixcr_airr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MiXCR_exportairr module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def clns_file = clns.size() > 1 ? clns.first() : clns // for sc input 2 clns are provided which create the same airr file. So, we take the first one.
    """
    # activate license
    if [ \${MIXCR_LICENSE:-"unset"} != "unset" ]; then
        echo "Initializing MIXCR_LICENSE env variable"
        export MI_LICENSE=\$MIXCR_LICENSE
    fi
    
    mixcr exportAirr \\
        ${clns_file} \\
        ${prefix}.airr.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcrexportairr: \$(mixcr --version |& sed '1!d ; s/mixcr //')
    END_VERSIONS
    """
}
