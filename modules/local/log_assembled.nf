process ASSEMBLED_LOGS {
    tag "logs"
    label 'process_low'

    conda "bioconda::pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'biocontainers/pandas:1.1.5' }"

    input:
    path('logs/*')
    path("p-sc-duplicates.csv")

    output:
    path "Table_sequences_assembled.tsv", emit: logs
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | grep -o "[0-9\\. ]\\+"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'pandas\').version)"'), emit: versions_pandas, topic: versions

    script:
    """
    log.py
    """
   
}