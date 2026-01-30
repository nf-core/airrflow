process PRESTO_PARSEHEADERS_PRIMERS {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_presto:0f8e73dc0555493d' :
        'biocontainers/presto:0.7.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    if (params.cprimer_position == "R1") {
        """
        ParseHeaders.py copy -s $reads -o ${reads.baseName}_reheader-pass.fastq -f $args --act first last -k C_PRIMER V_PRIMER

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( ParseHeaders.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    } else if (params.cprimer_position == "R2") {
        """
        ParseHeaders.py copy -s $reads -o ${reads.baseName}_reheader-pass.fastq -f $args --act first last -k V_PRIMER C_PRIMER

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            presto: \$( ParseHeaders.py --version | awk -F' '  '{print \$2}' )
        END_VERSIONS
        """
    }

}
