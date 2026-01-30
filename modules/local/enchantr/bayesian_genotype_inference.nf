def asString (args) {
    def s = ""
    def value = ""
    if (args.size()>0) {
        if (args[0] != 'none') {
            for (param in args.keySet().sort()){
                value = args[param].toString()
                if (!value.isNumber()) {
                    value = "'"+value+"'"
                }
                s = s + ",'"+param+"'="+value
            }
        }
    }
    return s
}

process BAYESIAN_GENOTYPE_INFERENCE {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    container "docker.io/immcantation/airrflow:genotyping"

    input:
    tuple val(meta), path(tabs) // meta, sequence tsv in AIRR format
    path reference_fasta
    path repertoires_samplesheet

    output:
    path "*_report/db_genotype", emit: reference // reference folder
    path("*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml", emit: versions


    script:
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = ""
    if (repertoires_samplesheet) {
        input = repertoires_samplesheet
    } else {
        input = tabs.join(',')
    }
    """
    Rscript -e "enchantr::enchantr_report('tigger_bayesian_genotype', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'species'='human', \\
                                        'genotypeby'='${params.genotypeby}', \\
                                        'single_clone_representative'='${params.single_clone_representative}', \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_bayesian_genotype_inference_command_log' ${args}))"

    cp -r enchantr ${meta.id}_bayesian_genotype_inference_report && rm -rf enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
