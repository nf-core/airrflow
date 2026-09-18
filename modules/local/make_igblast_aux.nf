process MAKE_IGBLAST_AUX {
    label 'process_low'

    // No conda directive: receptor_utils is pip-only, it is not on bioconda.
    // The sourcerer image carries receptor_utils and makeblastdb, which the
    // internal_data annotation database rebuild needs.
    container "docker.io/peresay/airrflow-sourcerer:sourcerer-bf1223f-5"

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
