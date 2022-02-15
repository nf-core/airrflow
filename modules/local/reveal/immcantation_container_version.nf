/*
 * Immcantation version
 */
process IMMCANTATION {
    label 'immcantation'
    label 'single_cpu'

    output:
    path "versions.yml", emit: versions

    script:
    //TODO: dev version is hardcoded now
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}_CONTAINER":
        immcantation: dev
    END_VERSIONS
    """
}
