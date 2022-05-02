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
    if ! command -v versions report &> /dev/null
    then
    cat <<-END_VERSIONS > versions.yml
    "${task.process}_CONTAINER":
        immcantation: none
    END_VERSIONS
    else
    echo "${task.process}_CONTAINER:" > versions.yml && \
    cat /Version.yaml | grep "^ " | grep -v "date:" | sed s/version/immcantation/g >> versions.yml
    fi
    """
}
