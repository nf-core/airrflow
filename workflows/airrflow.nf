/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Report files
ch_report_rmd       = Channel.fromPath(params.report_rmd, checkIfExists: true)
ch_report_css       = Channel.fromPath(params.report_css, checkIfExists: true)
ch_report_logo      = Channel.fromPath(params.report_logo, checkIfExists: true)
ch_report_logo_img  = Channel.fromPath(params.report_logo_img, checkIfExists: true)

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
include { REPERTOIRE_ANALYSIS_REPORTING } from '../subworkflows/local/repertoire_analysis_reporting'
include { SC_RAW_INPUT                  } from '../subworkflows/local/sc_raw_input'
include { MIXCR_FLOW                    } from '../subworkflows/local/mixcr_flow'
include { MIXCR_POSTANALYSIS            } from '../subworkflows/local/mixcr_postanalysis'
include { FASTQ_INPUT_CHECK             } from '../subworkflows/local/fastq_input_check'
include { RNASEQ_INPUT                  } from '../subworkflows/local/rnaseq_input'

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
        ch_input // channel: samplesheet read in from --input

    main:

        ch_versions = Channel.empty()
        ch_reassign_logs = Channel.empty()

        // Download or fetch databases
        DATABASES()

        if ( params.mode == "fastq" ) {

            // SC: Perform sequence assembly if input type is fastq from single-cell sequencing data (currently only 10XGenomics)
            if (params.library_generation_method == "sc_10x_genomics") {

                SC_RAW_INPUT(
                    ch_input
                )

                ch_fasta                                = SC_RAW_INPUT.out.fasta
                ch_versions                             = ch_versions.mix(SC_RAW_INPUT.out.versions)
                ch_cellranger_airr                      = SC_RAW_INPUT.out.airr
                ch_cellranger_out                       = SC_RAW_INPUT.out.outs

                ch_validated_samplesheet                = SC_RAW_INPUT.out.samplesheet.collect()

            ch_presto_filterseq_logs                = Channel.empty()
            ch_presto_maskprimers_logs              = Channel.empty()
            ch_presto_pairseq_logs                  = Channel.empty()
            ch_presto_clustersets_logs              = Channel.empty()
            ch_presto_buildconsensus_logs           = Channel.empty()
            ch_presto_postconsensus_pairseq_logs    = Channel.empty()
            ch_presto_assemblepairs_logs            = Channel.empty()
            ch_presto_collapseseq_logs              = Channel.empty()
            ch_presto_splitseq_logs                 = Channel.empty()
            ch_fastp_html                           = Channel.empty()
            ch_fastp_json                           = Channel.empty()
            ch_fastqc_postassembly_mqc              = Channel.empty()


        }  else if (params.library_generation_method == "trust4") {
            // Extract VDJ sequences from "general" RNA seq data using TRUST4

            RNASEQ_INPUT (
                ch_input,
                DATABASES.out.igblast.collect()
            )

            ch_fasta                                = RNASEQ_INPUT.out.fasta
            ch_versions                             = ch_versions.mix(RNASEQ_INPUT.out.versions)

            ch_validated_samplesheet                = RNASEQ_INPUT.out.samplesheet.collect()

            ch_presto_filterseq_logs                = Channel.empty()
            ch_presto_maskprimers_logs              = Channel.empty()
            ch_presto_pairseq_logs                  = Channel.empty()
            ch_presto_clustersets_logs              = Channel.empty()
            ch_presto_buildconsensus_logs           = Channel.empty()
            ch_presto_postconsensus_pairseq_logs    = Channel.empty()
            ch_presto_assemblepairs_logs            = Channel.empty()
            ch_presto_collapseseq_logs              = Channel.empty()
            ch_presto_splitseq_logs                 = Channel.empty()
            ch_fastp_html                           = RNASEQ_INPUT.out.fastp_reads_html
            ch_fastp_json                           = RNASEQ_INPUT.out.fastp_reads_json
            ch_fastqc_postassembly_mqc              = Channel.empty()
        }
        else if (params.library_generation_method == "mixcr") {

                if (!params.kit) {
                    error "Kit parameter is required for MiXCR analysis."
                }

                MIXCR_FLOW(ch_input)

                ch_fasta                    = MIXCR_FLOW.out.fasta
                ch_versions                 = ch_versions.mix(MIXCR_FLOW.out.versions)
                ch_mixcr_airr               = MIXCR_FLOW.out.airr
                ch_mixcr_clns               = MIXCR_FLOW.out.clns
                ch_mixcr_out                = MIXCR_FLOW.out.outs

                ch_validated_samplesheet = MIXCR_FLOW.out.samplesheet.collect()

                ch_presto_filterseq_logs             = Channel.empty()
                ch_presto_maskprimers_logs           = Channel.empty()
                ch_presto_pairseq_logs               = Channel.empty()
                ch_presto_clustersets_logs           = Channel.empty()
                ch_presto_buildconsensus_logs        = Channel.empty()
                ch_presto_postconsensus_pairseq_logs = Channel.empty()
                ch_presto_assemblepairs_logs         = Channel.empty()
                ch_presto_collapseseq_logs           = Channel.empty()
                ch_presto_splitseq_logs              = Channel.empty()
                ch_fastp_html                        = Channel.empty()
                ch_fastp_json                        = Channel.empty()
                ch_fastqc_postassembly_mqc           = Channel.empty()
            }
            else {
            // Perform sequence assembly if input type is fastq from bulk sequencing data
            SEQUENCE_ASSEMBLY(
                ch_input,
                DATABASES.out.igblast.collect()
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
            }

        } else if ( params.mode == "assembled" ) {

            ASSEMBLED_INPUT_CHECK (
                ch_input,
                params.miairr,
                params.collapseby,
                params.cloneby
            )
            ch_versions = ch_versions.mix( ASSEMBLED_INPUT_CHECK.out.versions )

            if (params.reassign) {
                CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                    ASSEMBLED_INPUT_CHECK.out.ch_tsv
                )
                ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta
                ch_versions = ch_versions.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.versions)
                ch_reassign_logs = ch_reassign_logs.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.logs)
            } else {
                ch_fasta_from_tsv = Channel.empty()
            }

            ch_fasta = ASSEMBLED_INPUT_CHECK.out.ch_fasta.mix(ch_fasta_from_tsv)
            ch_validated_samplesheet = ASSEMBLED_INPUT_CHECK.out.validated_input.collect()

            ch_presto_filterseq_logs             = Channel.empty()
            ch_presto_maskprimers_logs           = Channel.empty()
            ch_presto_pairseq_logs               = Channel.empty()
            ch_presto_clustersets_logs           = Channel.empty()
            ch_presto_buildconsensus_logs        = Channel.empty()
            ch_presto_postconsensus_pairseq_logs = Channel.empty()
            ch_presto_assemblepairs_logs         = Channel.empty()
            ch_presto_collapseseq_logs           = Channel.empty()
            ch_presto_splitseq_logs              = Channel.empty()
            ch_fastp_html                        = Channel.empty()
            ch_fastp_json                        = Channel.empty()
            ch_fastqc_postassembly_mqc           = Channel.empty()

        } else {
            error "Mode parameter value not valid."
        }
        // Perform V(D)J annotation and filtering
        VDJ_ANNOTATION(
            ch_fasta,
            ch_validated_samplesheet.collect(),
            DATABASES.out.igblast.collect(),
            DATABASES.out.reference_fasta.collect()
        )
        ch_versions = ch_versions.mix( VDJ_ANNOTATION.out.versions )

        // Split bulk and single cell repertoires
        ch_repertoire_by_processing = VDJ_ANNOTATION.out.repertoire
            .branch { it ->
                single: it[0].single_cell == 'true'
                bulk:   it[0].single_cell == 'false'
            }

        // Bulk: Assign germlines and filtering
        ch_repertoire_by_processing.bulk
            .dump(tag: 'bulk')

        BULK_QC_AND_FILTER(
            ch_repertoire_by_processing.bulk,
            VDJ_ANNOTATION.out.reference_fasta.collect()
        )
        ch_versions = ch_versions.mix( BULK_QC_AND_FILTER.out.versions )

        ch_bulk_filtered = BULK_QC_AND_FILTER.out.repertoires

        // Single cell: QC and filtering
        ch_repertoire_by_processing.single
            .dump(tag: 'single')

        SINGLE_CELL_QC_AND_FILTERING(
            ch_repertoire_by_processing.single
        )
        ch_versions = ch_versions.mix( SINGLE_CELL_QC_AND_FILTERING.out.versions )

        // Mixing bulk and single cell channels for clonal analysis
        ch_repertoires_for_clones = ch_bulk_filtered
                                        .mix(SINGLE_CELL_QC_AND_FILTERING.out.repertoires)
                                        .dump(tag: 'sc bulk mix')

        // Clonal analysis
        CLONAL_ANALYSIS(
            ch_repertoires_for_clones,
            VDJ_ANNOTATION.out.reference_fasta.collect(),
            ch_report_logo_img.collect().ifEmpty([])
        )
        ch_versions = ch_versions.mix( CLONAL_ANALYSIS.out.versions)

        if (!params.skip_report){
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
                ch_reassign_logs.collect().ifEmpty([]),
                VDJ_ANNOTATION.out.changeo_makedb_logs.collect().ifEmpty([]),
                VDJ_ANNOTATION.out.logs.collect().ifEmpty([]),
                BULK_QC_AND_FILTER.out.logs.collect().ifEmpty([]),
                SINGLE_CELL_QC_AND_FILTERING.out.logs.collect().ifEmpty([]),
                CLONAL_ANALYSIS.out.logs.collect().ifEmpty([]),
                CLONAL_ANALYSIS.out.repertoire,
                ch_input.collect(),
                ch_report_rmd.collect(),
                ch_report_css.collect(),
                ch_report_logo.collect(),
                ch_validated_samplesheet.collect()
            )
        }
        ch_versions = ch_versions.mix( REPERTOIRE_ANALYSIS_REPORTING.out.versions )

    // MiXCR postanalysis
    if (params.mixcr_postanalysis) {
        MIXCR_POSTANALYSIS ( ch_mixcr_clns )
    }


    // Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


        // MODULE: MultiQC

        if (!params.skip_multiqc) {
            summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

            ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
            ch_methods_description  = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

            ch_multiqc_files = Channel.empty()
            ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_json.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_postassembly_mqc.collect{it[1]}.ifEmpty([]))

            MULTIQC (
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_report_logo.toList(),
                [],
                []
            )
            multiqc_report = MULTIQC.out.report.toList()
        } else {
            multiqc_report = Channel.empty()
        }
    emit:
        multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
