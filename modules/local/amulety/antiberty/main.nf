process AMULETY_ANTIBERTY {
    tag "${meta.id}"
    label (params.use_gpu? 'process_gpu': 'process_medium')
    label 'immcantation'

    container "quay.io/biocontainers/mulled-v2-92ebbfc09fc136b8e201cb187cd9567ba335d439:459e6ebe51fb2818cb6de807f2c5fa99599b1214-0"
    publishDir "${params.outdir}/amulety/antiberty/", mode: 'copy'

    input:
    tuple val(meta), path(tsv) // meta, sequence tsv in AIRR format

    output:
    tuple val(meta), path("${tsv.baseName}_antiberty.tsv") , emit: embedding
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    amulety antiberty ${args} $tsv ${params.embedding_mode} ${tsv.baseName}_antiberty.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
    END_VERSIONS
    """
}
