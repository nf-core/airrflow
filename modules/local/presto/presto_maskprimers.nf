process PRESTO_MASKPRIMERS {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda "bioconda::presto=0.7.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c538a0e310c303233164cbe486b5e5f6bddcf18975a9b20ac2f590f151f03e62/data' :
        'biocontainers/presto:0.7.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)
    path(cprimers)
    path(vprimers)

    output:
    tuple val(meta), path("*_R1_primers-pass.fastq"), path("*_R2_primers-pass.fastq") , emit: reads
    path "*_command_log_R?.txt", emit: logs
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions


    script:
    def revpr = params.primer_revpr ? '--revpr' : ''
    if (params.cprimer_position == "R1") {
        def primer_start_R1 = (params.index_file | params.umi_position == 'R1') ? "--start ${params.umi_length + params.cprimer_start} --barcode" : "--start ${params.cprimer_start}"
        def primer_start_R2 = (params.umi_position == 'R2') ? "--start ${params.umi_length + params.vprimer_start} --barcode" : "--start ${params.vprimer_start}"
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $R1 -p ${cprimers} $primer_start_R1 $revpr --maxerror ${params.primer_r1_maxerror} --mode ${params.primer_r1_mask_mode} --outname ${meta.id}_R1 --log ${meta.id}_R1.log > ${meta.id}_command_log_R1.txt
        MaskPrimers.py score --nproc ${task.cpus} -s $R2 -p ${vprimers} $primer_start_R2 $revpr --maxerror ${params.primer_r2_maxerror} --mode ${params.primer_r2_mask_mode} --outname ${meta.id}_R2 --log ${meta.id}_R2.log > ${meta.id}_command_log_R2.txt
        ParseLog.py -l ${meta.id}_R1.log ${meta.id}_R2.log -f ID PRIMER ERROR

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else if (params.cprimer_position == "R2") {
        def primer_start_R1 = (params.index_file | params.umi_position == 'R1') ? "--start ${params.umi_length + params.vprimer_start} --barcode" : "--start ${params.vprimer_start}"
        def primer_start_R2 = (params.umi_position == 'R2') ? "--start ${params.umi_length + params.cprimer_start} --barcode" : "--start ${params.cprimer_start}"
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $R1 -p ${vprimers} $primer_start_R1 $revpr --maxerror ${params.primer_r1_maxerror} --mode ${params.primer_r1_mask_mode} --outname ${meta.id}_R1 --log ${meta.id}_R1.log > ${meta.id}_command_log_R1.txt
        MaskPrimers.py score --nproc ${task.cpus} -s $R2 -p ${cprimers} $primer_start_R2 $revpr --maxerror ${params.primer_r2_maxerror} --mode ${params.primer_r2_mask_mode} --outname ${meta.id}_R2 --log ${meta.id}_R2.log > ${meta.id}_command_log_R2.txt
        ParseLog.py -l "${meta.id}_R1.log" "${meta.id}_R2.log" -f ID PRIMER ERROR

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else {
        error "Error in determining cprimer position. Please choose R1 or R2."
    }

}
