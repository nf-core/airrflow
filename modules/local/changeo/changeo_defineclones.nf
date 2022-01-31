process CHANGEO_DEFINECLONES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::changeo=1.2.0 bioconda::igblast=1.17.1" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:f479475bceae84156e57e303cfe804ab5629d62b-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    val(threshold) // threshold file
    path(geno_fasta) // igblast fasta

    output:
    tuple val(meta), path("*clone-pass.tsv"), emit: tab // sequence tsv table in AIRR format
    path "*_command_log.txt" , emit: logs
    path("${geno_fasta}"), emit: fasta // genotype fasta

    script:
    if (params.set_cluster_threshold) {
        thr = params.cluster_threshold
    } else {
        thr = file(threshold).text
        thr = thr.trim()
    }
    """
    DefineClones.py -d $tab --act set --model ham --norm len --nproc $task.cpus --dist $thr --outname ${meta.id} --log ${meta.id}.log > "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}.log" -f id v_call j_call junction_length cloned filtered clones
    """
}
