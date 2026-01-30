process CHANGEO_CONVERTDB_FASTA {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'


    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/changeo_igblast_wget:dcfe290eb28df215' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*.fasta"), emit: fasta // sequence in fasta format
    path "versions.yml" , emit: versions
    path "*_command_log.txt" , emit: logs

    script:
    def args = task.ext.args ?: ''
    """
    ConvertDb.py fasta -d $tab $args > ${meta.id}_convertdb_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( ConvertDb.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
