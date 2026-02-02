// Import generic module functions
process RENAME_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}.${file.extension}"), emit: file
    path("*_command_log.txt"), emit: logs

    script:
    """
    mv ${file} ${meta.id}.${file.extension}
    echo "START> RenameFile" > ${meta.id}_rename_command_log.txt
    echo "FILE> ${file}" >> ${meta.id}_rename_command_log.txt
    echo "OUTPUT> ${meta.id}.${file.extension}" >> ${meta.id}_rename_command_log.txt
    if [[ "${file.extension}" == "fasta" ]]; then
        seq_count=\$(grep -c "^>" ${meta.id}.${file.extension})
        echo "RECORDS> \${seq_count}" >> ${meta.id}_rename_command_log.txt
        echo "PASS> \${seq_count}" >> ${meta.id}_rename_command_log.txt
    fi
    if [[ "${file.extension}" == "tsv" ]]; then
        seq_count=\$(( \$( grep -c '.' "${meta.id}.${file.extension}" ) - 1))
        echo "RECORDS> \${seq_count}" >> ${meta.id}_rename_command_log.txt
        echo "PASS> \${seq_count}" >> ${meta.id}_rename_command_log.txt
    fi
    """
}
