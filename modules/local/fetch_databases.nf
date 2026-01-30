process FETCH_DATABASES {
    tag "IMGT IGBLAST"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/changeo_igblast_wget:dcfe290eb28df215' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

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
    fetch_databases.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        IMGT download date: \$( echo \$(date "+%F") )
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
        changeo: \$( AssignGenes.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
