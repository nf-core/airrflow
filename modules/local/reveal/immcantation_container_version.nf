/*
 * Immcantation version
 */
process IMMCANTATION {
    label 'immcantation'
    label 'single_cpu'

    output:
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    if ! command -v versions report &> /dev/null
    then
        echo "immcantation: none" > ${software}.version.txt
    else
        versions report | head -n 1 > ${software}.version.txt
    fi
    """
}
