// Import generic module functions
include { saveFiles; getSoftwareName } from '../functions'

params.options = [:]

/*
 * Generate fasta from from AIRR rearrangement tsv file
 */
process FASTA_FROM_TSV {
    tag "${input_id}"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'fasta_files', publish_id:'') }

    conda (params.enable_conda ? { exit 1 "TODO: set up conda for FASTA_FROM_TSV." } : null)
    if (params.immcantation_container) {
         container params.immcantation_container
    } else if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         // container "https://depot.galaxyproject.org/singularity/mulled-v2-2665a8a48fa054ad1fcccf53e711669939b3eac1:09e1470e7d75ed23a083425eb01ce0418c9e8827-0"  // Singularity image
    } else {
        // container "quay.io/biocontainers/python:3.8.3"  // Docker image
    }

    input:
    tuple file(filename), val(subject_id), val(organism), val(collapseby_group), val(collapseby_size), val(cloneby_group), val(cloneby_size), val(filetype), val(input_id)
    
    output:
    tuple file("*.fasta"), val(subject_id), val(organism), val(collapseby_group), val(collapseby_size), val(cloneby_group), val(cloneby_size), val(filetype), val(input_id), emit: ch_fasta_from_tsv
    path "*.version.txt" , emit: version
       
    script:
    def software = getSoftwareName(task.process)
    """
    ConvertDb.py fasta -d "${filename}" --if sequence_id --sf sequence --mf cell_id consensus_count duplicate_count c_call c_cigar c_sequence_start c_sequence_end
    ConvertDb.py --version >  FASTA_FROM_TSV.version.txt
    """
}
