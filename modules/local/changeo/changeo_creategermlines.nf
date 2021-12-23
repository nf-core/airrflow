process CHANGEO_CREATEGERMLINES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' :
        'quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0' }"

    input:
    tuple val(meta), path(tab) // sequence tsv table in AIRR format
    path(geno_fasta) // igblast fasta
    path(imgt_base) // imgt db

    output:
    tuple val(meta), path("*germ-pass.tsv"), emit: tab
    path("*_command_log.txt"), emit: logs

    script:
    def software = getSoftwareName(task.process)
    """
    CreateGermlines.py -d ${tab} -g dmask --cloned -r ${geno_fasta} \\
    ${imgt_base}/${meta.species}/vdj/ \\
    --log ${meta.id}.log --outname ${meta.id} > ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}.log -f ID V_CALL D_CALL J_CALL
    """
}
