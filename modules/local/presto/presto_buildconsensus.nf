process PRESTO_BUILDCONSENSUS {
    tag "$meta.id"
    label "process_long_parallelized"
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::presto=0.7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.0--pyhdfd78af_0' :
        'quay.io/biocontainers/presto:0.7.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*_R1_consensus-pass.fastq"), path("*_R2_consensus-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*_R1.log")
    path("*_R2.log")
    path("*.tab"), emit: log_tab
    path "versions.yml" , emit: versions

    script:
    def barcode_field = params.cluster_sets ? "CLUSTER" : "BARCODE"
    """
    BuildConsensus.py -s $R1 --bf ${barcode_field} --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror $params.buildconsensus_maxerror --maxgap $params.buildconsensus_maxgap --outname ${meta.id}_R1 --log ${meta.id}_R1.log > "${meta.id}_command_log.txt"
    BuildConsensus.py -s $R2 --bf ${barcode_field} --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror $params.buildconsensus_maxerror --maxgap $params.buildconsensus_maxgap --outname ${meta.id}_R2 --log ${meta.id}_R2.log >> "${meta.id}_command_log.txt"
    ParseLog.py -l "${meta.id}_R1.log" "${meta.id}_R2.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( BuildConsensus.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
