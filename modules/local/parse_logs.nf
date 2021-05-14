// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PARSE_LOGS {
    tag "logs"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"logs") }

    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

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
    path('igblast/*') //CHANGEO_PARSEDB_SELECT logs
    path('define_clones/*') //CHANGEO_DEFINECLONES logs
    path('create_germlines/*') //CHANGEO_CREATEGERMLINES logs
    path('metadata.tsv') //METADATA

    output:
    path "Table_sequences_process.tsv", emit: logs

    script:
    '''
    log_parsing.py
    '''
}
