process PRESTO_ASSEMBLEPAIRS_SEQUENTIAL {
    tag "$meta.id"
    label 'process_long_parallelized'
    label 'immcantation'

    conda "bioconda::presto=0.7.8 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0 conda-forge::biopython=1.85"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/43dc39c7cc72fea3c2ad6615c19919587b8a413598b4044e93119e34b7a4d5eb/data':
        'community.wave.seqera.io/library/igblast_presto_biopython_wget:688a053e4779d278' }"

    input:
    tuple val(meta), path(R1), path(R2) // reads in fastq format
    path(igblast) // igblast references

    output:
    tuple val(meta), path("*_assemble-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*_table.tab")
    tuple val("${task.process}"), val('presto'), eval('AssemblePairs.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    AssemblePairs.py sequential -1 $R2 -2 $R1 --nproc ${task.cpus} \\
        -r "${igblast}/fasta/${meta.species}_${meta.locus.toLowerCase()}_v.fasta" \\
        $args \\
        --outname ${meta.id} --log ${meta.id}.log > ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}.log $args2

    """
}
