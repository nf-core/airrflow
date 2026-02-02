process PRESTO_ASSEMBLEPAIRS_SEQUENTIAL {
    tag "$meta.id"
    label 'process_long_parallelized'
    label 'immcantation'

    conda "bioconda::presto=0.7.6 bioconda::igblast=1.22.0 conda-forge::wget=1.21.4 conda-forge::biopython=1.85"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/633f1ff1426784e70bde4cd52ea1e46f690fd97876e4fbedd8dfe01977f94e5a/data':
        'community.wave.seqera.io/library/igblast_presto_biopython_wget:82da8ebc3ab7e416' }"

    input:
    tuple val(meta), path(R1), path(R2) // reads in fastq format
    path(igblast) // igblast references

    output:
    tuple val(meta), path("*_assemble-pass.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*_table.tab")
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    AssemblePairs.py sequential -1 $R2 -2 $R1 --nproc ${task.cpus} \\
        -r "${igblast}/fasta/imgt_${meta.species}_${meta.locus.toLowerCase()}_v.fasta" \\
        $args \\
        --outname ${meta.id} --log ${meta.id}.log > ${meta.id}_command_log.txt
    ParseLog.py -l ${meta.id}.log $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( AssemblePairs.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
