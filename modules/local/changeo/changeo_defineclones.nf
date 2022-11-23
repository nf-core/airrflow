process CHANGEO_DEFINECLONES {
    tag "$meta.id"
    label 'process_medium'
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::changeo=1.2.0 bioconda::igblast=1.17.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    val(threshold) // threshold file

    output:
    tuple val(meta), path("*clone-pass.tsv"), emit: tab // sequence tsv table in AIRR format
    path "*_command_log.txt" , emit: logs
    path "versions.yml" , emit: versions

    script:
    if (params.set_cluster_threshold) {
        thr = params.cluster_threshold
    } else {
        thr = file(threshold).text
        thr = thr.trim()
    }
    """
    DefineClones.py -d $tab --act set --model ham --norm len --nproc $task.cpus --dist $thr --outname ${meta.id} --log ${meta.id}.log > "${meta.id}_defineclones_command_log.txt"
    ParseLog.py -l "${meta.id}.log" -f id v_call j_call junction_length cloned filtered clones

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( DefineClones.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
