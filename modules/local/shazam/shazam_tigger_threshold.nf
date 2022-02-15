process SHAZAM_TIGGER_THRESHOLD {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::r-base=4.1.2 bioconda::r-shazam=1.1.0 conda-forge::r-tigger=1.0.0" : null)              // Conda package
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-db80433cc6df75bc43a5fd7bfa7529a7df8cfe15:f0e1329252bbc0f36a8656cfa655cf205da30e5b-0' :
        'quay.io/biocontainers/mulled-v2-db80433cc6df75bc43a5fd7bfa7529a7df8cfe15:f0e1329252bbc0f36a8656cfa655cf205da30e5b-0' }"

    input:
    tuple val(meta), path(tab) // tsv tab in AIRR format
    path(imgt_base) // igblast fasta

    output:
    tuple val(meta), path("*genotyped.tab"), emit: tab
    path("*threshold.txt"), emit: threshold
    path("versions.yml") , emit: versions
    path("*genotype.fasta"), emit: fasta
    path("*genotype.pdf")
    path("*Hamming_distance_threshold.pdf")

    script:
    def args = task.ext.args ?: ''
    def locus = meta.locus
    def species = meta.species
    if ( locus.equals('IG')){
        """
        TIgGER-shazam.R $tab ${locus.toLowerCase()} $params.threshold_method \
            "${imgt_base}/${species}/vdj/imgt_${species}_IGHV.fasta" \
            "${imgt_base}/${species}/vdj/imgt_${species}_IGKV.fasta" \
            "${imgt_base}/${species}/vdj/imgt_${species}_IGLV.fasta"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tigger: \$(Rscript -e "library(tigger); cat(paste(packageVersion('tigger'), collapse='.'))")
            shazam: \$(Rscript -e "library(shazam); cat(paste(packageVersion('shazam'), collapse='.'))")
            R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
        END_VERSIONS
        """
    } else if ( locus.equals('TR')){
        """
        TIgGER-shazam.R $tab ${locus.toLowerCase()} $params.threshold_method \\
        "${imgt_base}/${species}/vdj/imgt_${species}_TRAV.fasta" \\
        "${imgt_base}/${species}/vdj/imgt_${species}_TRBV.fasta" \\
        "${imgt_base}/${species}/vdj/imgt_${species}_TRDV.fasta" \\
        "${imgt_base}/${species}/vdj/imgt_${species}_TRGV.fasta"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tigger: \$(Rscript -e "library(tigger); cat(paste(packageVersion('tigger'), collapse='.'))")
            shazam: \$(Rscript -e "library(shazam); cat(paste(packageVersion('shazam'), collapse='.'))")
            R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
        END_VERSIONS
        """
    } else {
        """
        echo "Locus not supported, or mixture of locus provided. Choose from: 'IG','TR'"
        echo $locus
        exit 1
        """
    }
}
