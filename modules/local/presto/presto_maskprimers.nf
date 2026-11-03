process PRESTO_MASKPRIMERS {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)
    path(cprimers)
    path(vprimers)
    val primer_revpr
    val cprimer_position
    val index_file
    val umi_position
    val umi_length
    val cprimer_start
    val vprimer_start
    val primer_r1_maxerror
    val primer_r1_mask_mode
    val primer_r2_maxerror
    val primer_r2_mask_mode

    output:
    tuple val(meta), path("*_R1_primers-pass.fastq"), path("*_R2_primers-pass.fastq") , emit: reads
    path "*_command_log_R?.txt", emit: logs
    path "*.tab", emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('MaskPrimers.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions


    script:
    def revpr = primer_revpr ? '--revpr' : ''
    if (cprimer_position == "R1") {
        def primer_start_R1 = (index_file | umi_position == 'R1') ? "--start ${umi_length + cprimer_start} --barcode" : "--start ${cprimer_start}"
        def primer_start_R2 = (umi_position == 'R2') ? "--start ${umi_length + vprimer_start} --barcode" : "--start ${vprimer_start}"
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $R1 -p ${cprimers} $primer_start_R1 $revpr --maxerror ${primer_r1_maxerror} --mode ${primer_r1_mask_mode} --outname ${meta.id}_R1 --log ${meta.id}_R1.log > ${meta.id}_command_log_R1.txt
        MaskPrimers.py score --nproc ${task.cpus} -s $R2 -p ${vprimers} $primer_start_R2 $revpr --maxerror ${primer_r2_maxerror} --mode ${primer_r2_mask_mode} --outname ${meta.id}_R2 --log ${meta.id}_R2.log > ${meta.id}_command_log_R2.txt
        ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log -f ID PRIMER ERROR

        """
    } else if (cprimer_position == "R2") {
        def primer_start_R1 = (index_file | umi_position == 'R1') ? "--start ${umi_length + vprimer_start} --barcode" : "--start ${vprimer_start}"
        def primer_start_R2 = (umi_position == 'R2') ? "--start ${umi_length + cprimer_start} --barcode" : "--start ${cprimer_start}"
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $R1 -p ${vprimers} $primer_start_R1 $revpr --maxerror ${primer_r1_maxerror} --mode ${primer_r1_mask_mode} --outname ${meta.id}_R1 --log ${meta.id}_R1.log > ${meta.id}_command_log_R1.txt
        MaskPrimers.py score --nproc ${task.cpus} -s $R2 -p ${cprimers} $primer_start_R2 $revpr --maxerror ${primer_r2_maxerror} --mode ${primer_r2_mask_mode} --outname ${meta.id}_R2 --log ${meta.id}_R2.log > ${meta.id}_command_log_R2.txt
        ParseLog.py -l "${meta.id}_R1.log" "${meta.id}_R2.log" -f ID PRIMER ERROR

        """
    } else {
        error "Error in determining cprimer position. Please choose R1 or R2."
    }

}
