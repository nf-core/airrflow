process FETCH_DATABASES {
    tag "$database_type IGBLAST"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "docker.io/peresay/airrflow-sourcerer:sourcerer-bf1223f-5"

    input:
    val(database_type)

    output:
    path("igblast_base"), emit: igblast
    path("reference_base"), emit: reference_fasta
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('changeo'), eval('AssignGenes.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_changeo, topic: versions
    path("igblast_base/database/human_ig_v.ndb"), emit: igblast_human_ig_v
    path("igblast_base/database/human_ig_d.ndb"), emit: igblast_human_ig_d
    path("igblast_base/database/human_ig_j.ndb"), emit: igblast_human_ig_j
    path("igblast_base/database/human_tr_v.ndb"), emit: igblast_human_tr_v
    path("igblast_base/database/human_tr_d.ndb"), emit: igblast_human_tr_d
    path("igblast_base/database/human_tr_j.ndb"), emit: igblast_human_tr_j

    script:
    """
    fetch_databases.sh -d ${database_type}

    """
}
