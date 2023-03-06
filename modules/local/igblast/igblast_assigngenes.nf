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
    tuple val(meta), path("*db-pass.tsv"), emit: tab
    path "versions.yml" , emit: versions
    path("*_command_log.txt"), emit: logs //process logs
    path("*_makedb_command_log.txt"), emit: makedb_log

    script:
    def args = task.ext.args ?: ''
    """
    igblastn \
        -germline_db_V igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_v \
        -germline_db_D igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_d \
        -germline_db_J igblast_base/database/imgt_${meta.species}_${meta.locus.toLowerCase()}_j \
        -auxiliary_data igblast_base/optional_file/${meta.species}_gl.aux \
        -organism ${meta.species} \
        $args \
        -query $reads \
        -out ${meta.id}_db-pass.tsv

    echo "START> AssignGenes" > ${meta.id}_changeo_assigngenes_command_log.txt
    echo "COMMAND> igblast" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "VERSION> \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "FILE> ${reads}" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "ORGANISM> ${meta.species}" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "LOCI> ${meta.locus.toLowerCase()}" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "NPROC> ${task.cpus}\n" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "PROGRESS> ...Done \n" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "PASS> \$(tail -n +2 ${meta.id}_db-pass.tsv | wc -l )" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "OUTPUT> ${meta.id}_igblast.fmt7" >> ${meta.id}_changeo_assigngenes_command_log.txt
    echo "END> AssignGenes" >> ${meta.id}_changeo_assigngenes_command_log.txt

    echo "START> MakeDB" > ${meta.id}_makedb_command_log.txt
    echo "COMMAND> igblast" >> ${meta.id}_makedb_command_log.txt
    echo "ALIGNER_FILE> ${meta.id}_igblast.fmt7" >> ${meta.id}_makedb_command_log.txt
    echo "SEQ_FILE> ${reads}" >> ${meta.id}_makedb_command_log.txt
    echo "ASIS_ID> False" >> ${meta.id}_makedb_command_log.txt
    echo "ASIS_CALLS> False" >> ${meta.id}_makedb_command_log.txt
    echo "VALIDATE> strict" >> ${meta.id}_makedb_command_log.txt
    echo "EXTENDED> True" >> ${meta.id}_makedb_command_log.txt
    echo "INFER_JUNCTION> False\n" >> ${meta.id}_makedb_command_log.txt
    echo "PROGRESS> ...\n" >> ${meta.id}_makedb_command_log.txt
    echo "PROGRESS> ... Done\n" >> ${meta.id}_makedb_command_log.txt
    echo "OUTPUT> ${meta.id}_db-pass.tsv" >> ${meta.id}_makedb_command_log.txt
    echo "PASS> \$(tail -n +2 ${meta.id}_db-pass.tsv | wc -l )" >> ${meta.id}_makedb_command_log.txt
    echo "FAIL> 0" >> ${meta.id}_makedb_command_log.txt
    echo "END> MakeDB" >> ${meta.id}_makedb_command_log.txt

    echo "\"${task.process}\":" > versions.yml
    echo "   igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )" >> versions.yml
    """
}
