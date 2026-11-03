process PRESTO_MASKPRIMERS_POSTASSEMBLY {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(cprimers)
    path(vprimers)
    val primer_revpr
    val cprimer_position
    val cprimer_start
    val primer_r1_maxerror
    val primer_r1_mask_mode
    val vprimer_start
    val primer_r2_maxerror
    val primer_r2_mask_mode

    output:
    tuple val(meta), path("*REV_primers-pass.fastq") , emit: reads
    path "*command_log.txt", emit: logs
    path "*.tab", emit: log_tab
    tuple val("${task.process}"), val('presto'), eval('MaskPrimers.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def revpr = primer_revpr ? '--revpr' : ''
    if (cprimer_position == "R1") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${cprimers} --start ${cprimer_start} --maxerror ${primer_r1_maxerror} \
            --mode ${primer_r1_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > ${meta.id}_command_log.txt
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${vprimers} --start ${vprimer_start} --maxerror ${primer_r2_maxerror} \
            --mode ${primer_r2_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> ${meta.id}_command_log.txt
        ParseLog.py -l ${meta.id}-FWD.log ${meta.id}-REV.log -f ID PRIMER ERROR

        """
    } else if (cprimer_position == "R2") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${vprimers} --start ${vprimer_start} --maxerror ${primer_r1_maxerror} \
            --mode ${primer_r1_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > ${meta.id}_command_log.txt
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${cprimers} --start ${cprimer_start} --maxerror ${primer_r2_maxerror} \
            --mode ${primer_r2_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> ${meta.id}_command_log.txt
        ParseLog.py -l ${meta.id}-FWD.log ${meta.id}-REV.log -f ID PRIMER ERROR

        """
    } else {
        error "Error in determining cprimer positon."
    }

}
