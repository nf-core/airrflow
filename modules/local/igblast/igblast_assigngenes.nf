process IGBLAST_ASSIGNGENES {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.0 bioconda::igblast=1.19.0 conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:00534555924705cdf2f7ac48b4b8b4083527ca58-1' :
        'quay.io/biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:00534555924705cdf2f7ac48b4b8b4083527ca58-1' }"

    input:
    tuple val(meta), path(reads) // reads in fasta format
    path(igblast) // igblast fasta

    output:
    path("*igblast.fmt7"), emit: blast
    tuple val(meta), path("$reads"), emit: fasta
    path "versions.yml" , emit: versions
    path("*_command_log.txt"), emit: logs //process logs

    script:
    def args = task.ext.args ?: ''
    """
    igblastn \
        -germline_db_V igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_v\
        -germline_db_D igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_d \
        -germline_db_J igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_j \
        -auxiliary_data igblast_base/optional_file/${meta.species}_gl.aux \
        -organism ${meta.species} \
        $args \
        -query $reads \
        -out ${meta.id}_igblast.fmt7

    echo "START> AssignGenes" > ${meta.id}_changeo_assigngenes_command_log.txt
    echo "COMMAND> igblast" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "FILE> ${reads}" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "PASS> \$(tail -n 1 ${meta.id}_igblast.fmt7 | grep -o "[0-9]\\+" )" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "END> AssignGenes" >> ${meta.id}_changeo_assigngenes_command_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
    END_VERSIONS
    """
}
