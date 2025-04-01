process AMULETY_TRANSLATE {
    tag "${meta.id}"
    label "process_medium"
    label 'immcantation'

    container "docker.io/immcantation/airrflow:4.3.0"
    publishDir "${params.outdir}/translations/", mode: "copy"

    input:
    tuple val(meta), path(tsv) // meta, sequence tsv in AIRR format
    path(igblast) // igblast references

    output:
    path "*_translated.tsv" , emit: translated

    script:
    """
    amulety translate-igblast $tsv . igblast
    """
}
