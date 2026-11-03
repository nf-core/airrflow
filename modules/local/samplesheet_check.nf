process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::pandas=1.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'biocontainers/pandas:1.1.5' }"

    input:
    path samplesheet

    output:
    path '*.tsv', emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | grep -o "[0-9\\. ]\\+"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'pandas\').version)"'), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/airrflow/bin/
    def args = task.ext.args ?: ''
    """
    check_samplesheet.py $samplesheet $args
    cp $samplesheet samplesheet.valid.tsv

    """
}
