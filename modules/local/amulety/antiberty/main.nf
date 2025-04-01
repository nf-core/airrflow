process AMULETY_ANTIBERTY{
    tag "${meta.id}"
    label (params.use_gpu? 'process_gpu': 'process_medium')
    label 'immcantation'

    container "docker.io/immcantation/airrflow:4.3.0"
    publishDir "${params.outdir}/amulety/antiberty/", mode: 'copy'

    input:
    tuple val(meta), path(tsv) // meta, sequence tsv in AIRR format

    output:
    path "${tsv.baseName}_antiberty.tsv" , emit: embedding

    script:
    def args = task.ext.args ?: ''
    """
    amulety antiberty ${args} $tsv ${params.embedding_mode} ${tsv.baseName}_antiberty.tsv
    """
}
