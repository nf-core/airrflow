// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
def options    = initOptions(params.options)

process CHANGEO_ASSIGNGENES {
    tag "${input_id}"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'assign_genes', publish_id:'') }

    conda (params.enable_conda ? "bioconda::changeo=1.0.2 bioconda::igblast=1.15.0" : null)              // Conda package
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        container "quay.io/biocontainers/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"                        // Docker image
    }

    input:
    tuple file(repertoire_fasta), val(subject_id), val(organism), val(collapseby_group), val(collapseby_size), val(cloneby_group), val(cloneby_size), val(filetype), val(input_id)
    val(igblast_db )
    val(imgt_db)
    val(igblastn)

    output:
    tuple file("*.fmt7"), val(subject_id), val(organism), val(collapseby_group), val(collapseby_size), val(cloneby_group), val(cloneby_size), val(filetype ), val(input_id), emit: ch_fmt7
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    AssignGenes.py igblast -s "${repertoire_fasta}" \
      -b $igblast_db --organism ${organism} --loci ig \
      --format blast --nproc $task.cpus --exec $igblastn
    AssignGenes.py --version > ASSIGN_GENES.version.txt
    """
}
