#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/airrflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/airrflow
    Website: https://nf-co.re/airrflow
    Slack  : https://nfcore.slack.com/channels/airrflow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AIRRFLOW                } from './workflows/airrflow'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_airrflow_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_airrflow_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_airrflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER DECLARATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params {

    // Input parameters
    input: String? = null
    mode: String = "fastq"
    miairr: String = "$projectDir/assets/reveal/mapping_MiAIRR_BioSample_v1.3.1.tsv"
    index_file: Boolean = false
    save_intermediate_repertoires: Boolean = true

    // ----------------------------
    // sequencing protocol options
    // ----------------------------
    library_generation_method: String? = null
    race_linker: String? = null

    // Primer and UMI inputs
    cprimers: String? = null
    vprimers: String? = null
    vprimer_start: Integer = 0
    cprimer_start: Integer = 0
    cprimer_position: String = 'R1'
    primer_revpr: Boolean = false

    // UMI and primer handling
    umi_position: String = 'R1'
    umi_length: Integer = -1
    umi_start: Integer = 0

    // trimming options
    trim_fastq: Boolean = true
    adapter_fasta: String? = null
    clip_r1: Integer = 0
    clip_r2: Integer = 0
    three_prime_clip_r1: Integer = 0
    three_prime_clip_r2: Integer = 0
    trim_nextseq: Boolean = false
    save_trimmed: Boolean = false

    // --------------------------
    // sequence assembly options
    // --------------------------
    // Filter sequences
    filterseq_q: Integer = 20

    // Mask primers
    primer_r1_maxerror: Float = 0.2
    primer_r2_maxerror: Float = 0.2
    primer_r1_mask_mode: String = 'cut'
    primer_r2_mask_mode: String = 'cut'
    maskprimers_align_race: Boolean = false
    maskprimers_align: Boolean = false
    maskprimers_extract: Boolean = false
    primer_r1_extract_len: Integer = 0
    primer_r2_extract_len: Integer = 0
    primer_maxlen: Integer = 50

    // Build consensus
    primer_consensus: Float = 0.6
    buildconsensus_maxerror: Float = 0.1
    buildconsensus_maxgap: Float = 0.5
    cluster_sets: Boolean = true

    // Assemble pairs
    assemblepairs_sequential: Boolean = false

    // internal cregion
    align_cregion: Boolean = false
    internal_cregion_sequences: String? = null
    cregion_maxlen: Integer = 100
    cregion_maxerror: Float = 0.3
    cregion_mask_mode: String = 'tag'

    // -----------------------
    // vdj annotation options
    // -----------------------
    productive_only: Boolean = true
    reassign: Boolean = true
    reference_igblast: String = params.pipelines_testdata_base_path + 'database-cache/igblast_base.zip'
    reference_fasta: String = params.pipelines_testdata_base_path + 'database-cache/imgtdb_base.zip'
    fetch_germlines: String = 'none'
    skip_alignment_filter: Boolean = false

    // -----------------------
    // bulk filtering options
    // -----------------------
    remove_chimeric: Boolean = false
    detect_contamination: Boolean = false
    collapseby: String = 'sample_id'

    // -----------------------
    // clonal analysis options
    // -----------------------
    cloneby: String = 'subject_id'
    crossby: String = 'subject_id'
    singlecell: String = 'single_cell'
    clonal_threshold = 'auto'
    skip_all_clones_report: Boolean = false
    skip_report_threshold: Boolean = false
    skip_clonal_analysis: Boolean = false

    // -----------------------
    // tree lineage options
    // -----------------------
    lineage_tree_builder: String = 'raxml'
    lineage_tree_exec: String = '/usr/local/bin/raxml-ng'
    lineage_trees: Boolean = false

    // -----------------------
    // genotyping options
    // -----------------------
    genotyping: Boolean = false
    genotypeby: String = 'subject_id'
    novel_allele_inference: Boolean = true
    genotyping_clonal_threshold: Float = 0.2
    single_clone_representative: Boolean = true

    // -----------------------
    // translate embed options
    // -----------------------
    translate: Boolean = false
    embeddings: String? = null
    embedding_chain: String = 'H'
    use_gpu: Boolean = false

    // -----------------------
    // reporting options
    // -----------------------
    skip_report: Boolean = false
    report_rmd: String = "$projectDir/assets/repertoire_comparison.Rmd"
    report_css: String = "$projectDir/assets/nf-core_style.css"
    report_logo: String = "$projectDir/assets/nf-core-airrflow_logo_light.png"
    report_logo_img: String = "$projectDir/assets/nf-core-airrflow_logo_reports.png"

    // -----------------------
    // Single cell raw input options
    // -----------------------
    reference_10x: String? = null

    // -----------------------
    // raw RNA seq input options
    // -----------------------
    trust4_cell_barcode_read: String? = null
    trust4_read_format: String? = null
    trust4_umi_read: String? = null
    trust4_barcode_whitelist: String? = null


    // -----------------------
    // generic nf-core options
    // -----------------------

    // References
    igenomes_base: String = 's3://ngi-igenomes/igenomes/'
    igenomes_ignore: Boolean = true

    // MultiQC options
    skip_multiqc: Boolean = false
    multiqc_config: String? = null
    multiqc_title: String? = null
    multiqc_logo: String? = null
    max_multiqc_email_size: String = '25.MB'
    multiqc_methods_description: String? = null

    // Boilerplate options
    outdir: String? = null
    email: String? = null
    email_on_fail: String? = null
    plaintext_email: Boolean = false
    monochrome_logs: Boolean = false
    help: Boolean = false
    help_full: Boolean = false
    show_hidden: Boolean = false
    version: Boolean = false

    // Schema validation default options
    validate_params: Boolean = true
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_AIRRFLOW {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? file( params.multiqc_config, checkIfExists: true ) : []
    ch_multiqc_logo          = params.multiqc_logo   ? file( params.multiqc_logo, checkIfExists: true ) : []
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // Report files
    ch_report_rmd       = channel.fromPath(params.report_rmd, checkIfExists: true)
    ch_report_css       = channel.fromPath(params.report_css, checkIfExists: true)
    ch_report_logo      = channel.fromPath(params.report_logo, checkIfExists: true)
    ch_report_logo_img  = channel.fromPath(params.report_logo_img, checkIfExists: true)

    //
    // WORKFLOW: Run pipeline
    //
    AIRRFLOW (
        samplesheet,
        params.mode,
        params.library_generation_method,
        params.miairr,
        params.collapseby,
        params.cloneby,
        params.reassign,
        params.genotyping,
        params.skip_clonal_analysis,
        params.translate,
        params.embeddings,
        params.skip_report,
        params.outdir,
        params.skip_multiqc,
        params.multiqc_methods_description,
        ch_report_rmd,
        ch_report_css,
        ch_report_logo,
        ch_report_logo_img,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        params.fetch_germlines,
        params.reference_igblast,
        params.reference_fasta,
        params.vprimers,
        params.race_linker,
        params.cprimers,
        params.umi_length,
        params.reference_10x,
        params.index_file,
        params.trust4_barcode_whitelist,
        params.trust4_cell_barcode_read,
        params.trust4_umi_read,
        params.trust4_read_format,
        params.skip_alignment_filter,
        params.productive_only,
        params.remove_chimeric,
        params.detect_contamination,
        params.genotypeby,
        params.novel_allele_inference,
        params.single_clone_representative,
        params.genotyping_clonal_threshold,
        params.clonal_threshold,
        params.skip_report_threshold,
        params.skip_all_clones_report,
        params.lineage_trees,
        params.embedding_chain,
        params.adapter_fasta,
        params.maskprimers_extract,
        params.internal_cregion_sequences,
        params.maskprimers_align_race,
        params.umi_position,
        params.umi_start,
        params.save_trimmed,
        params.maskprimers_align,
        params.cprimer_position,
        params.primer_maxlen,
        params.primer_r1_maxerror,
        params.primer_r1_mask_mode,
        params.primer_r2_maxerror,
        params.primer_r2_mask_mode,
        params.cprimer_start,
        params.vprimer_start,
        params.primer_revpr,
        params.primer_r2_extract_len,
        params.primer_r1_extract_len,
        params.cluster_sets,
        params.assemblepairs_sequential,
        params.align_cregion,
        params.cregion_maxlen,
        params.cregion_maxerror,
        params.cregion_mask_mode,
        params.crossby,
        params.singlecell,
        params.lineage_tree_builder,
        params.lineage_tree_exec,
        params.filterseq_q,
        params.buildconsensus_maxerror,
        params.buildconsensus_maxgap,
        params.primer_consensus
    )
    emit:
    multiqc_report = AIRRFLOW.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_AIRRFLOW (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_AIRRFLOW.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
