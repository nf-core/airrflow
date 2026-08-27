/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR } from '../modules/local/changeo/changeo_convertdb_fasta'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { DATABASES                     } from '../subworkflows/local/databases'
include { SEQUENCE_ASSEMBLY             } from '../subworkflows/local/sequence_assembly'
include { ASSEMBLED_INPUT_CHECK         } from '../subworkflows/local/assembled_input_check'
include { VDJ_ANNOTATION                } from '../subworkflows/local/vdj_annotation'
include { BULK_QC_AND_FILTER            } from '../subworkflows/local/bulk_qc_and_filter'
include { SINGLE_CELL_QC_AND_FILTERING  } from '../subworkflows/local/single_cell_qc_and_filtering'
include { CLONAL_ANALYSIS               } from '../subworkflows/local/clonal_analysis'
include { NOVEL_ALLELES_AND_GENOTYPING   } from '../subworkflows/local/novel_alleles_and_genotyping'
include { REPERTOIRE_ANALYSIS_REPORTING } from '../subworkflows/local/repertoire_analysis_reporting'
include { SC_RAW_INPUT                  } from '../subworkflows/local/sc_raw_input'
include { RNASEQ_INPUT                  } from '../subworkflows/local/rnaseq_input'
include { TRANSLATE_EMBED              } from '../subworkflows/local/translate_embed'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_airrflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AIRRFLOW {

    take:
        ch_input
        mode
        library_generation_method
        miairr
        collapseby
        cloneby
        reassign
        genotyping
        skip_clonal_analysis
        translate
        embeddings
        skip_report
        outdir
        skip_multiqc
        multiqc_methods_description
        ch_report_rmd
        ch_report_css
        ch_report_logo
        ch_report_logo_img
        ch_multiqc_config
        ch_multiqc_custom_config
        ch_multiqc_logo
        fetch_germlines
        reference_igblast
        reference_fasta
        vprimers
        race_linker
        cprimers
        umi_length
        reference_10x
        index_file
        trust4_barcode_whitelist
        trust4_cell_barcode_read
        trust4_umi_read
        trust4_read_format
        skip_alignment_filter
        productive_only
        remove_chimeric
        detect_contamination
        genotypeby
        novel_allele_inference
        single_clone_representative
        genotyping_clonal_threshold
        clonal_threshold
        skip_report_threshold
        skip_all_clones_report
        lineage_trees
        embedding_chain
        adapter_fasta
        maskprimers_extract
        internal_cregion_sequences
        maskprimers_align_race
        umi_position
        umi_start
        save_trimmed
        maskprimers_align
        cprimer_position
        primer_maxlen
        primer_r1_maxerror
        primer_r1_mask_mode
        primer_r2_maxerror
        primer_r2_mask_mode
        cprimer_start
        vprimer_start
        primer_revpr
        primer_r2_extract_len
        primer_r1_extract_len
        cluster_sets
        assemblepairs_sequential
        align_cregion
        cregion_maxlen
        cregion_maxerror
        cregion_mask_mode
        crossby
        singlecell
        lineage_tree_builder
        lineage_tree_exec
        filterseq_q
        buildconsensus_maxerror
        buildconsensus_maxgap
        primer_consensus

    main:

        ch_versions = channel.empty()
        ch_reassign_logs = channel.empty()
        ch_input_check_logs = channel.empty()

        // Download or fetch databases
        DATABASES(
            fetch_germlines,
            reference_igblast,
            reference_fasta
        )

        if ( mode == "fastq" ) {

            // SC:Perform sequence assembly if input type is fastq from single-cell sequencing data (currently only 10XGenomics)
            if (library_generation_method == "sc_10x_genomics") {

                SC_RAW_INPUT(
                    ch_input,
                    vprimers,
                    race_linker,
                    cprimers,
                    umi_length,
                    reference_10x,
                    library_generation_method,
                    collapseby,
                    cloneby,
                    index_file
                )

                ch_fasta                                = SC_RAW_INPUT.out.fasta
                ch_versions                             = ch_versions.mix(SC_RAW_INPUT.out.versions)
                ch_cellranger_airr                      = SC_RAW_INPUT.out.airr
                ch_cellranger_out                       = SC_RAW_INPUT.out.outs
                ch_validated_samplesheet                = SC_RAW_INPUT.out.samplesheet.collect()
                ch_presto_filterseq_logs                = channel.empty()
                ch_presto_maskprimers_logs              = channel.empty()
                ch_presto_pairseq_logs                  = channel.empty()
                ch_presto_clustersets_logs              = channel.empty()
                ch_presto_buildconsensus_logs           = channel.empty()
                ch_presto_postconsensus_pairseq_logs    = channel.empty()
                ch_presto_assemblepairs_logs            = channel.empty()
                ch_presto_collapseseq_logs              = channel.empty()
                ch_presto_splitseq_logs                 = channel.empty()
                ch_fastp_html                           = channel.empty()
                ch_fastp_json                           = channel.empty()
                ch_fastqc_postassembly_mqc              = channel.empty()
                ch_tsv_files                            = channel.empty()

            }  else if (library_generation_method == "trust4") {
                // Extract VDJ sequences from "general" RNA seq data using TRUST4

                RNASEQ_INPUT (
                    ch_input,
                    DATABASES.out.igblast.collect(),
                    vprimers,
                    race_linker,
                    cprimers,
                    umi_length,
                    reference_10x,
                    trust4_barcode_whitelist,
                    trust4_cell_barcode_read,
                    trust4_umi_read,
                    trust4_read_format,
                    library_generation_method,
                    collapseby,
                    cloneby,
                    index_file
                )

                ch_fasta                                = RNASEQ_INPUT.out.fasta
                ch_versions                             = ch_versions.mix(RNASEQ_INPUT.out.versions)
                ch_validated_samplesheet                = RNASEQ_INPUT.out.samplesheet.collect()

                ch_presto_filterseq_logs                = channel.empty()
                ch_presto_maskprimers_logs              = channel.empty()
                ch_presto_pairseq_logs                  = channel.empty()
                ch_presto_clustersets_logs              = channel.empty()
                ch_presto_buildconsensus_logs           = channel.empty()
                ch_presto_postconsensus_pairseq_logs    = channel.empty()
                ch_presto_assemblepairs_logs            = channel.empty()
                ch_presto_collapseseq_logs              = channel.empty()
                ch_presto_splitseq_logs                 = channel.empty()
                ch_fastp_html                           = RNASEQ_INPUT.out.fastp_reads_html
                ch_fastp_json                           = RNASEQ_INPUT.out.fastp_reads_json
                ch_fastqc_postassembly_mqc              = channel.empty()
                ch_tsv_files                            = channel.empty()
            } else {
                // Perform sequence assembly if input type is fastq from bulk sequencing data
                SEQUENCE_ASSEMBLY(
                    ch_input,
                    DATABASES.out.igblast.collect(),
                    library_generation_method,
                    adapter_fasta,
                    maskprimers_extract,
                    vprimers,
                    cprimers,
                    race_linker,
                    umi_length,
                    internal_cregion_sequences,
                    maskprimers_align_race,
                    index_file,
                    umi_position,
                    umi_start,
                    collapseby,
                    cloneby,
                    save_trimmed,
                    maskprimers_align,
                    cprimer_position,
                    primer_maxlen,
                    primer_r1_maxerror,
                    primer_r1_mask_mode,
                    primer_r2_maxerror,
                    primer_r2_mask_mode,
                    cprimer_start,
                    vprimer_start,
                    primer_revpr,
                    primer_r2_extract_len,
                    primer_r1_extract_len,
                    cluster_sets,
                    assemblepairs_sequential,
                    align_cregion,
                    cregion_maxlen,
                    cregion_maxerror,
                    cregion_mask_mode,
                    filterseq_q,
                    buildconsensus_maxerror,
                    buildconsensus_maxgap,
                    primer_consensus
                )

                ch_fasta                                = SEQUENCE_ASSEMBLY.out.fasta
                ch_versions                             = ch_versions.mix(SEQUENCE_ASSEMBLY.out.versions)
                ch_fastp_html                           = SEQUENCE_ASSEMBLY.out.fastp_reads_html
                ch_fastp_json                           = SEQUENCE_ASSEMBLY.out.fastp_reads_json
                ch_fastqc_postassembly_mqc              = SEQUENCE_ASSEMBLY.out.fastqc_postassembly
                ch_validated_samplesheet                = SEQUENCE_ASSEMBLY.out.samplesheet.collect()
                ch_presto_filterseq_logs                = SEQUENCE_ASSEMBLY.out.presto_filterseq_logs.ifEmpty([])
                ch_presto_maskprimers_logs              = SEQUENCE_ASSEMBLY.out.presto_maskprimers_logs.ifEmpty([])
                ch_presto_pairseq_logs                  = SEQUENCE_ASSEMBLY.out.presto_pairseq_logs.ifEmpty([])
                ch_presto_clustersets_logs              = SEQUENCE_ASSEMBLY.out.presto_clustersets_logs.ifEmpty([])
                ch_presto_buildconsensus_logs           = SEQUENCE_ASSEMBLY.out.presto_buildconsensus_logs.ifEmpty([])
                ch_presto_postconsensus_pairseq_logs    = SEQUENCE_ASSEMBLY.out.presto_postconsensus_pairseq_logs.ifEmpty([])
                ch_presto_assemblepairs_logs            = SEQUENCE_ASSEMBLY.out.presto_assemblepairs_logs.ifEmpty([])
                ch_presto_collapseseq_logs              = SEQUENCE_ASSEMBLY.out.presto_collapseseq_logs.ifEmpty([])
                ch_presto_splitseq_logs                 = SEQUENCE_ASSEMBLY.out.presto_splitseq_logs.ifEmpty([])
                ch_tsv_files                            = channel.empty()
            }

        } else if ( mode == "assembled" ) {

            ASSEMBLED_INPUT_CHECK (
                ch_input,
                miairr,
                collapseby,
                cloneby,
                reassign
            )
            ch_input_check_logs = ASSEMBLED_INPUT_CHECK.out.logs

            if (reassign) {
                CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                    ASSEMBLED_INPUT_CHECK.out.ch_tsv
                )
                ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta
                ch_reassign_logs = ch_reassign_logs.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.logs)
                ch_tsv_files = channel.empty()
            } else {
                ch_fasta_from_tsv = channel.empty()
                ch_tsv_files = ASSEMBLED_INPUT_CHECK.out.ch_tsv
            }

            ch_fasta = ASSEMBLED_INPUT_CHECK.out.ch_fasta.mix(ch_fasta_from_tsv)
            ch_validated_samplesheet = ASSEMBLED_INPUT_CHECK.out.validated_input.collect()

            ch_presto_filterseq_logs             = channel.empty()
            ch_presto_maskprimers_logs           = channel.empty()
            ch_presto_pairseq_logs               = channel.empty()
            ch_presto_clustersets_logs           = channel.empty()
            ch_presto_buildconsensus_logs        = channel.empty()
            ch_presto_postconsensus_pairseq_logs = channel.empty()
            ch_presto_assemblepairs_logs         = channel.empty()
            ch_presto_collapseseq_logs           = channel.empty()
            ch_presto_splitseq_logs              = channel.empty()
            ch_fastp_html                        = channel.empty()
            ch_fastp_json                        = channel.empty()
            ch_fastqc_postassembly_mqc           = channel.empty()

        } else {
            error "Mode parameter value not valid."
        }

        // Perform V(D)J annotation and filtering
        VDJ_ANNOTATION(
            ch_fasta,
            ch_tsv_files,
            ch_validated_samplesheet.collect(),
            DATABASES.out.igblast.collect(),
            DATABASES.out.reference_fasta.collect(),
            skip_alignment_filter,
            productive_only
        )

        // Split bulk and single cell repertoires
        ch_repertoire_by_processing = VDJ_ANNOTATION.out.repertoire
            .branch { it ->
                single: it[0].single_cell == 'true'
                bulk:   it[0].single_cell == 'false'
            }

        // Bulk: Assign germlines and filtering
        ch_repertoire_by_processing.bulk

        BULK_QC_AND_FILTER(
            ch_repertoire_by_processing.bulk,
            VDJ_ANNOTATION.out.reference_fasta.collect(),
            remove_chimeric,
            detect_contamination,
            collapseby
        )

        ch_bulk_filtered = BULK_QC_AND_FILTER.out.repertoires

        // Single cell: QC and filtering
        ch_repertoire_by_processing.single

        SINGLE_CELL_QC_AND_FILTERING(
            ch_repertoire_by_processing.single
        )

        // Mixing bulk and single cell channels after filtering
        ch_repertoires_after_qc = ch_bulk_filtered
                                        .mix(SINGLE_CELL_QC_AND_FILTERING.out.repertoires)

        // Novel alleles and genotype inference
        if (genotyping) {
            NOVEL_ALLELES_AND_GENOTYPING(
                ch_repertoires_after_qc,
                VDJ_ANNOTATION.out.reference_fasta.collect(),
                ch_validated_samplesheet.collect(),
                ch_report_logo_img.collect().ifEmpty([]),
                genotypeby,
                novel_allele_inference,
                single_clone_representative,
                genotyping_clonal_threshold,
                cloneby,
                singlecell
            )
            ch_repertoire_reference = NOVEL_ALLELES_AND_GENOTYPING.out.repertoire_reference

        } else {
            ch_repertoire_reference = ch_repertoires_after_qc.combine(VDJ_ANNOTATION.out.reference_fasta)
        }
        ch_repertoire_reference.dump(tag: 'ch_repertoire_reference_forcloning')

        // Clonal analysis
        if (!skip_clonal_analysis) {
            CLONAL_ANALYSIS(
                ch_repertoire_reference,
                ch_report_logo_img.collect().ifEmpty([]),
                clonal_threshold,
                skip_report_threshold,
                cloneby,
                skip_all_clones_report,
                lineage_trees,
                genotypeby,
                crossby,
                singlecell,
                lineage_tree_builder,
                lineage_tree_exec
            )
        }

        // Translation and embedding
        if (translate || embeddings) {
            TRANSLATE_EMBED(
                ch_repertoires_after_qc,
                DATABASES.out.igblast.collect(),
                embeddings,
                embedding_chain
            )
        }

        if (!skip_report){

            REPERTOIRE_ANALYSIS_REPORTING(
                ch_presto_filterseq_logs.collect().ifEmpty([]),
                ch_presto_maskprimers_logs.collect().ifEmpty([]),
                ch_presto_pairseq_logs.collect().ifEmpty([]),
                ch_presto_clustersets_logs.collect().ifEmpty([]),
                ch_presto_buildconsensus_logs.collect().ifEmpty([]),
                ch_presto_postconsensus_pairseq_logs.collect().ifEmpty([]),
                ch_presto_assemblepairs_logs.collect().ifEmpty([]),
                ch_presto_collapseseq_logs.collect().ifEmpty([]),
                ch_presto_splitseq_logs.collect().ifEmpty([]),
                ch_input_check_logs.collect().ifEmpty([]),
                ch_reassign_logs.collect().ifEmpty([]),
                VDJ_ANNOTATION.out.changeo_makedb_logs.collect().ifEmpty([]),
                VDJ_ANNOTATION.out.logs.collect().ifEmpty([]),
                BULK_QC_AND_FILTER.out.logs.collect().ifEmpty([]),
                SINGLE_CELL_QC_AND_FILTERING.out.logs.collect().ifEmpty([]),
                ch_input.collect(),
                ch_report_rmd.collect(),
                ch_report_css.collect(),
                ch_report_logo.collect(),
                ch_validated_samplesheet.collect(),
                mode,
                library_generation_method,
                umi_length,
                cluster_sets
            )
        }


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            def escaped_version = version.toString().replace('\\', '\\\\').replace('"', '\\"')
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: \"${escaped_version}\"" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'airrflow_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

        // MODULE: MultiQC

        if (!skip_multiqc) {
            summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

            ch_multiqc_custom_methods_description = multiqc_methods_description ? file(multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
            ch_methods_description  = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

            ch_multiqc_files = channel.empty()
            ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_json.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_postassembly_mqc.collect{it[1]}.ifEmpty([]))

            ch_multiqc_files_collected = ch_multiqc_files
                .collect()
                .ifEmpty([])

            // Build MultiQC input tuple
            def buildMultiqcInputTuple = { id, files ->
                [
                    [id: id],
                    files,
                    [ch_multiqc_config, ch_multiqc_custom_config].findAll { cfg -> cfg },
                    ch_multiqc_logo,
                    [],
                    []
                ]
            }

            // Merge all multiqc input channels
            ch_multiqc_input = ch_multiqc_files_collected
                .map { files ->
                    buildMultiqcInputTuple.call('multiqc_report', files)
                }

            ch_multiqc_input.dump(tag: 'ch_multiqc_input_before_multiqc')

            MULTIQC (
                ch_multiqc_input
            )
            multiqc_report = MULTIQC.out.report.toList()
        } else {
            multiqc_report = channel.empty()
        }
    emit:
        multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
