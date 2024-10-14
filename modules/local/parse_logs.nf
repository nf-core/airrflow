// Import generic module functions
process PARSE_LOGS {
    tag "logs"
    label 'process_low'

    conda "bioconda::pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'biocontainers/pandas:1.1.5' }"

    input:
    path('filter_by_sequence_quality/*') //PRESTO_FILTERSEQ logs
    path('mask_primers/*') //PRESTO_MASKPRIMERS logs
    path('pair_sequences/*') //PRESTO_PAIRSEQ logs
    path('cluster_sets/*') //PRESTO_CLUTSERSETS logs
    path('build_consensus/*') //PRESTO_BUILDCONSENSUS logs
    path('repair_mates/*') //PRESTO_POSTCONSESUS_PAIRSEQ logs
    path('assemble_pairs/*') //PRESTO_ASSEMBLEPAIRS logs
    path('deduplicates/*') //PRESTO_COLLAPSESEQ logs
    path('filter_representative_2/*') //PRESTO_SPLITSEQ logs
    path('igblast/*') //CHANGEO_MAKEDB logs
    path('metadata.tsv') //METADATA

    output:
    path "Table_sequences_process.tsv", emit: logs
    path "Table*.tsv", emit:tables
    path "versions.yml" , emit: versions

    script:
    if (params.umi_length == 0) {
        """
        log_parsing_no-umi.py

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$( echo \$(python --version | grep -o "[0-9\\. ]\\+") )
            pandas: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)"))
        END_VERSIONS
        """
    } else {
        def clustersets = params.cluster_sets? "--cluster_sets":""
        """
        log_parsing.py $clustersets

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$( echo \$(python --version | grep -o "[0-9\\. ]\\+") )
            pandas: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)"))
        END_VERSIONS
        """
    }
}
