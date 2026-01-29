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

process REASSIGN_ALLELES {
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
    val segments // which segments to reassign alleles to 
    //TODO: did we want to handle all segments at once? Then this val channel would not be needed.
    // *After novel alleles we just need to change the V, it's a time waste to go over all segments.
    //TODO: Check if we need the outputby parameter. Right now this is the same as the cloneby parameter.
    output:
    path("*/*/db_genotype"), emit: reference // reference folder
    path("*/*_reassigned.tsv"), emit: repertoires // reassigned repertoire
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    path "versions.yml", emit: versions


    script:
    def args = task.ext.args ? asString(task.ext.args) : ''
    def segs = segments.join(",")
    def input = ""
    if (repertoires_samplesheet) {
        input = repertoires_samplesheet
    } else {
        input = tabs.join(',')
    }
    """
    Rscript -e "enchantr::enchantr_report('reassign_alleles', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'species'='auto', \\
                                        'outputby'='${params.cloneby}', \\
                                        'segments'='${segs}', \\
                                        'force'=FALSE, \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_reassign_alleles_command_log' ${args}))"

    cp -r enchantr ${meta.id}_reassign_alleles_report && rm -rf enchantr

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  enchantr: ',packageVersion('enchantr'),'\n'))" >> versions.yml
    """
}
