process MAKE_IGBLAST_AUX {
    label 'process_low'

    // blast provides makeblastdb, which the internal_data V database rebuild needs.
    conda "bioconda::blast=2.17.0 conda-forge::python=3.11 pip:receptor-utils==0.0.67"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/2540941dbd0f691b5f5cb49cb37118e1efefd538897954e28d05be42e0d17fbe/data' :
        'community.wave.seqera.io/library/blast_python_pip_receptor-utils:9470ae020c3cb9d0' }"

    input:
    path(igblast_dir, stageAs: "input_igblast_base")
    path(reference_fasta_dir)

    output:
    path("igblast_base"), emit: igblast
    tuple val("${task.process}"), val('receptor_utils'), eval('python3 -c \'from importlib.metadata import version; print(version("receptor_utils"))\''), emit: versions_receptor_utils, topic: versions

    script:
    """
    make_igblast_aux.py -i input_igblast_base -r "${reference_fasta_dir}" -o igblast_base
    """
}
