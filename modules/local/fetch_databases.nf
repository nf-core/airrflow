process FETCH_DATABASES {
    tag "$database_type IGBLAST"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::igblast=1.22.0 bioconda::blast=2.17.0 conda-forge::python=3.11 pip:immvoke==0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e6e85f6cc5e8a5938c697dd3db88ce61fdeeb4ea7431a34937610e6918e459f/data' :
        'community.wave.seqera.io/library/blast_igblast_python_pip_immvoke:0293d51b5a9b5b88' }"

    input:
    val(database_type)

    output:
    path("igblast_base"), emit: igblast
    path("reference_base"), emit: reference_fasta
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+"'), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('immvoke'), eval('immvoke --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_immvoke, topic: versions

    script:
    """
    fetch_databases.sh -d ${database_type}
    """
}
