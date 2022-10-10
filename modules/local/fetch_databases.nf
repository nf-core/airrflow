process FETCH_DATABASES {
    tag "IMGT IGBLAST"
    label 'process_low'
    label 'immcantation'

    conda (params.enable_conda ? "bioconda::changeo=1.2.0 bioconda::igblast=1.19.0 conda-forge::wget=1.20.3" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e7f88c6f7da46a5407f261ca406c050d5bd12dea-0' :
        'quay.io/biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:e7f88c6f7da46a5407f261ca406c050d5bd12dea-0' }"

    output:
    path("igblast_base"), emit: igblast
    path("imgtdb_base"), emit: imgt
    path "versions.yml" , emit: versions
    path("igblast_base/database/imgt_human_ig_v.ndb"), emit: igblast_human_ig_v
    path("igblast_base/database/imgt_human_ig_d.ndb"), emit: igblast_human_ig_d
    path("igblast_base/database/imgt_human_ig_j.ndb"), emit: igblast_human_ig_j
    path("igblast_base/database/imgt_human_tr_v.ndb"), emit: igblast_human_tr_v
    path("igblast_base/database/imgt_human_tr_d.ndb"), emit: igblast_human_tr_d
    path("igblast_base/database/imgt_human_tr_j.ndb"), emit: igblast_human_tr_j

    script:
    """
    fetch_databases.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMGT download date: \$( echo \$(date "+%F") )
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
