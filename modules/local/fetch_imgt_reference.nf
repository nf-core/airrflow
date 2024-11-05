process FETCH_IMGT_REFERENCE {
    tag "IMGT reference"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.0 bioconda::igblast=1.22.0 conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' :
        'biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' }"

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: reference_fasta
    path "versions.yml" , emit: versions
    path("igblast_base/database/imgt_human_ig_v.ndb"), emit: igblast_human_ig_v
    path("igblast_base/database/imgt_human_ig_d.ndb"), emit: igblast_human_ig_d
    path("igblast_base/database/imgt_human_ig_j.ndb"), emit: igblast_human_ig_j
    path("igblast_base/database/imgt_human_tr_v.ndb"), emit: igblast_human_tr_v
    path("igblast_base/database/imgt_human_tr_d.ndb"), emit: igblast_human_tr_d
    path("igblast_base/database/imgt_human_tr_j.ndb"), emit: igblast_human_tr_j

    script:
    """
    fetch_imgt_reference.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMGT download date: \$( echo \$(date "+%F") )
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
