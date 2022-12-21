process PRESTO_MASKPRIMERS_POSTASSEMBLY {
    tag "$meta.id"
    label "process_high"
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::presto=0.7.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(cprimers)
    path(vprimers)

    output:
    tuple val(meta), path("*REV_primers-pass.fastq") , emit: reads
    path "*command_log.txt", emit: logs
    path "*.log"
    path "*.tab", emit: log_tab
    path "versions.yml" , emit: versions

    script:
    def revpr = params.primer_revpr ? '--revpr' : ''
    if (params.cprimer_position == "R1") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${cprimers} --start ${params.cprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > "${meta.id}_command_log.txt"
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${vprimers} --start ${params.vprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> "${meta.id}_command_log.txt"
        ParseLog.py -l "${meta.id}-FWD.log" "${meta.id}-REV.log" -f ID PRIMER ERROR

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else if (params.cprimer_position == "R2") {
        """
        MaskPrimers.py score --nproc ${task.cpus} -s $reads -p ${vprimers} --start ${params.cprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-FWD \
            --log ${meta.id}-FWD.log > "${meta.id}_command_log.txt"
        MaskPrimers.py score --nproc ${task.cpus} -s ${meta.id}-FWD_primers-pass.fastq -p ${cprimers} --start ${params.vprimer_start} --maxerror ${params.primer_maxerror} \
            --mode ${params.primer_mask_mode} --outname ${meta.id}-REV $revpr \
            --log ${meta.id}-REV.log >> "${meta.id}_command_log.txt"
        ParseLog.py -l "${meta.id}-FWD.log" "${meta.id}-REV.log" -f ID PRIMER ERROR

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( MaskPrimers.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else {
        exit 1, "Error in determining cprimer positon."
    }

}
