process MIXCR_IND_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    secret 'MIXCR_LICENSE'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0':
        'ghcr.io/milaboratory/mixcr/mixcr:4.6.0' }"

    input:
    tuple val(meta), path(mixcr_ind_json)
    val(diversity_plottype)
    val(statistical_method)
    val(p_adjust_method)

    output:
    tuple val(meta), path('*')  , emit: outs

    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MiXCR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # activate license
    if [ \${MIXCR_LICENSE:-"unset"} != "unset" ]; then
        echo "Initializing MIXCR_LICENSE env variable"
        export MI_LICENSE=\$MIXCR_LICENSE
    fi

    # individual diversity plots
    mixcr exportPlots diversity \\
        --plot-type ${diversity_plottype} \\
        --method ${statistical_method} \\
        --p-adjust-method ${p_adjust_method} \\
        ${mixcr_ind_json} \\
        ${prefix}.diversity.pdf

    # individual cdr3 plots
    mixcr exportPlots cdr3metrics \\
        ${mixcr_ind_json} \\
        --method ${statistical_method} \\
        --p-adjust-method ${p_adjust_method} \\
        ${prefix}.cdr3.pdf

    # V usage
    mixcr exportPlots vUsage \\
        ${mixcr_ind_json} \\
        ${prefix}.V_usage_heatmap.pdf
    
    mixcr exportPlots vUsage \\
        --bar-plot \\
        ${mixcr_ind_json} \\
        ${prefix}.V_usage_barplot.pdf
    
    mixcr exportPlots vUsage \\
        --bar-plot-by-sample \\
        ${mixcr_ind_json} \\
        ${prefix}.V_usage_barplot_by_sample.pdf

    # J usage
    mixcr exportPlots jUsage \\
        ${mixcr_ind_json} \\
        ${prefix}.J_usage_heatmap.pdf

    mixcr exportPlots jUsage \\
        --bar-plot \\
        ${mixcr_ind_json} \\
        ${prefix}.J_usage_barplot.pdf
    
    mixcr exportPlots jUsage \\
        --bar-plot-by-sample \\
        ${mixcr_ind_json} \\
        ${prefix}.J_usage_barplot_by_sample.pdf

    # VJ usage
    mixcr exportPlots vjUsage \\
        ${mixcr_ind_json} \\
        ${prefix}.VJ_usage_heatmap.pdf

    # Isotype usage
    mixcr exportPlots isotypeUsage \\
        ${mixcr_ind_json} \\
        ${prefix}.isotype_usage_heatmap.pdf

    mixcr exportPlots isotypeUsage \\
        --bar-plot \\
        ${mixcr_ind_json} \\
        ${prefix}.isotype_usage_barplot.pdf
    
    mixcr exportPlots isotypeUsage \\
        --bar-plot-by-sample \\
        ${mixcr_ind_json} \\
        ${prefix}.isotype_usage_barplot_by_sample.pdf


    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr --version |& sed '1!d ; s/mixcr //')
    END_VERSIONS
    """
}
