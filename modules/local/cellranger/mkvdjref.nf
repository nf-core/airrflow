process CELLRANGER_MKVDJREF {
    tag "${meta.id}"
    label 'process_high'

    container "nf-core/cellranger:7.1.0"

    input:
    path fasta
    // path gtf // maybe this can be used to decide whether --seqs or --genes should be used. If its possible to set gtf as an optional parameter.
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKVDJREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    cellranger \\
        mkvdjref \\
        --genome=$reference_name \\
        --seqs=$fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
