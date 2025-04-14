process AMULETY_TRANSLATE {
    tag "${meta.id}"
    label "process_medium"
    label 'immcantation'

    container "biocontainers/mulled-v2-92ebbfc09fc136b8e201cb187cd9567ba335d439:459e6ebe51fb2818cb6de807f2c5fa99599b1214-0"
    publishDir "${params.outdir}/translations/", mode: "copy"

    input:
    tuple val(meta), path(tsv) // meta, sequence tsv in AIRR format
    path(reference_igblast) // igblast references

    output:
    tuple val(meta), path("*_translated.tsv") , emit: repertoire_translated
    path "versions.yml" , emit: versions

    script:
    """
    export IGDATA=${reference_igblast}
    amulety translate-igblast $tsv . $reference_igblast

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
    END_VERSIONS
    """
}
