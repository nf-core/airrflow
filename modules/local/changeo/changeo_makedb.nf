include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process CHANGEO_MAKEDB {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
    }

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect()
    path(imgt_base)

    output:
    tuple val(meta), path("*db-pass.tsv"), emit: tab //sequence table in AIRR format
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def software = getSoftwareName(task.process)
    if (params.loci == 'ig'){
        """
        MakeDb.py igblast -i $igblast -s $reads -r \\
        ${imgt_base}/${params.species}/vdj/imgt_${params.species}_IG*.fasta \\
        --regions default --format airr --outname "${meta.id}" > "${meta.id}_${task.process}_command_log.txt"
        """
    } else if (params.loci == 'tr') {
        """
        MakeDb.py igblast -i $igblast -s $reads -r \\
        "${imgt_base}/${params.species}/vdj/imgt_${params.species}_TR*.fasta" \\
        --regions default --format airr --outname "${meta.id}" > "${meta.id}_${task.process}_command_log.txt"
        """
    }
}
