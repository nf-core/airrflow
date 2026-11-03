process PRESTO_PARSEHEADERS_PRIMERS {
    tag "$meta.id"
    label "process_low"
    label 'immcantation'

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val cprimer_position

    output:
    tuple val(meta), path("*_reheader-pass.fastq"), emit: reads
    tuple val("${task.process}"), val('presto'), eval('ParseHeaders.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args ?: ''
    if (cprimer_position == "R1") {
        """
        ParseHeaders.py copy -s $reads -o ${reads.baseName}_reheader-pass.fastq -f $args --act first last -k C_PRIMER V_PRIMER

        """
    } else if (cprimer_position == "R2") {
        """
        ParseHeaders.py copy -s $reads -o ${reads.baseName}_reheader-pass.fastq -f $args --act first last -k V_PRIMER C_PRIMER

        """
    }

}
