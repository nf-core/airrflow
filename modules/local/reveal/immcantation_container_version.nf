/*
 * Immcantation version
 */
process IMMCANTATION {
    label 'immcantation'
    label 'single_cpu'

    output:
    path "*.version.txt", emit: version

    script:
    """
    if ! command -v versions report &> /dev/null
    then
        echo "immcantation: none" > immcantation_container.version.txt
    else
        versions report | head -n 1 > immcantation_container.version.txt
    fi
    """
}
