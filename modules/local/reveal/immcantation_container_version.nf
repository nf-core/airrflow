/*
 * Immcantation version
 */
process IMMCANTATION {
    label 'immcantation'
    label 'single_cpu'

    output:
    path "versions.yml", emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}_CONTAINER":
        \$(echo \$(versions report) | head -n 1 )
    END_VERSIONS
    """
}
