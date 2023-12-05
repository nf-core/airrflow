/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAirrflow.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input) {
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
} else {
    error "Please provide input file containing the sample metadata with the '--input' option."
}


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
ch_report_rmd = Channel.fromPath(params.report_rmd, checkIfExists: true)
ch_report_css = Channel.fromPath(params.report_css, checkIfExists: true)
ch_report_logo = Channel.fromPath(params.report_logo, checkIfExists: true)
ch_report_logo_img = Channel.fromPath(params.report_logo_img, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR } from '../modules/local/changeo/changeo_convertdb_fasta'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SEQUENCE_ASSEMBLY             } from '../subworkflows/local/sequence_assembly'
include { ASSEMBLED_INPUT_CHECK         } from '../subworkflows/local/assembled_input_check'
include { VDJ_ANNOTATION                } from '../subworkflows/local/vdj_annotation'
include { BULK_QC_AND_FILTER            } from '../subworkflows/local/bulk_qc_and_filter'
include { SINGLE_CELL_QC_AND_FILTERING  } from '../subworkflows/local/single_cell_qc_and_filtering'
include { CLONAL_ANALYSIS               } from '../subworkflows/local/clonal_analysis'
include { REPERTOIRE_ANALYSIS_REPORTING } from '../subworkflows/local/repertoire_analysis_reporting'
include { SC_RAW_INPUT                  } from '../subworkflows/local/sc_raw_input'
include { FASTQ_INPUT_CHECK             } from '../subworkflows/local/fastq_input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow AIRRFLOW {

    ch_versions = Channel.empty()
    ch_reassign_logs = Channel.empty()

    if ( params.mode == "fastq" ) {

        // check the samplesheet to distinguish between bulk and single cell
        FASTQ_INPUT_CHECK(ch_input)
        ch_versions = ch_versions.mix(FASTQ_INPUT_CHECK.out.versions)

        ch_reads = FASTQ_INPUT_CHECK.out.reads

        // Split bulk and single cell repertoires
        ch_reads
            .branch { meta, reads ->
                single: meta.single_cell == 'true'
                bulk:   meta.single_cell == 'false'
            }
            .set { ch_reads_split }

        // SC: Perform downstream analysis
        SC_RAW_INPUT(
            ch_reads_split.single
        )


        // // Perform sequence assembly if input type is fastq
        // // TODO make this part run from ch_reads_split.bulk! -> other input, FASTQ_INPUT_CHECK is not needed then anymore

        // SEQUENCE_ASSEMBLY( ch_input )

        // ch_fasta                    = SEQUENCE_ASSEMBLY.out.fasta
        // ch_versions                 = ch_versions.mix(SEQUENCE_ASSEMBLY.out.versions)
        // ch_fastp_html               = SEQUENCE_ASSEMBLY.out.fastp_reads_html
        // ch_fastp_json               = SEQUENCE_ASSEMBLY.out.fastp_reads_json
        // ch_fastqc_postassembly_mqc  = SEQUENCE_ASSEMBLY.out.fastqc_postassembly
        // ch_validated_samplesheet    = SEQUENCE_ASSEMBLY.out.samplesheet.collect()

        // ch_presto_filterseq_logs        = SEQUENCE_ASSEMBLY.out.presto_filterseq_logs
        // ch_presto_maskprimers_logs      = SEQUENCE_ASSEMBLY.out.presto_maskprimers_logs
        // ch_presto_pairseq_logs          = SEQUENCE_ASSEMBLY.out.presto_pairseq_logs
        // ch_presto_clustersets_logs      = SEQUENCE_ASSEMBLY.out.presto_clustersets_logs
        // ch_presto_buildconsensus_logs   = SEQUENCE_ASSEMBLY.out.presto_buildconsensus_logs
        // ch_presto_postconsensus_pairseq_logs = SEQUENCE_ASSEMBLY.out.presto_postconsensus_pairseq_logs
        // ch_presto_assemblepairs_logs    = SEQUENCE_ASSEMBLY.out.presto_assemblepairs_logs
        // ch_presto_collapseseq_logs      = SEQUENCE_ASSEMBLY.out.presto_collapseseq_logs
        // ch_presto_splitseq_logs         = SEQUENCE_ASSEMBLY.out.presto_splitseq_logs

    } else if ( params.mode == "assembled" ) {

        ASSEMBLED_INPUT_CHECK (
            ch_input,
            params.miairr,
            params.collapseby,
            params.cloneby
        )
        ch_versions = ch_versions.mix( ASSEMBLED_INPUT_CHECK.out.versions.ifEmpty([]) )

        if (params.reassign) {
            CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                ASSEMBLED_INPUT_CHECK.out.ch_tsv
            )
            ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta
            ch_versions = ch_versions.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.versions.ifEmpty([]))
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
    // // Perform V(D)J annotation and filtering
    // VDJ_ANNOTATION(
    //     ch_fasta,
    //     ch_validated_samplesheet.collect()
    // )
    // ch_versions = ch_versions.mix( VDJ_ANNOTATION.out.versions.ifEmpty([]))

    // // Split bulk and single cell repertoires
    // ch_repertoire_by_processing = VDJ_ANNOTATION.out.repertoire
    //     .branch { it ->
    //         single: it[0].single_cell == 'true'
    //         bulk:   it[0].single_cell == 'false'
    //     }

    // // Bulk: Assign germlines and filtering
    // ch_repertoire_by_processing.bulk
    //     .dump(tag: 'bulk')

    // BULK_QC_AND_FILTER(
    //     ch_repertoire_by_processing.bulk,
    //     VDJ_ANNOTATION.out.imgt.collect()
    // )
    // ch_versions = ch_versions.mix( BULK_QC_AND_FILTER.out.versions.ifEmpty([]))

    // ch_bulk_filtered = BULK_QC_AND_FILTER.out.repertoires

    // // Single cell: QC and filtering
    // ch_repertoire_by_processing.single
    //     .dump(tag: 'single')

    // SINGLE_CELL_QC_AND_FILTERING(
    //     ch_repertoire_by_processing.single
    // )
    // ch_versions = ch_versions.mix( SINGLE_CELL_QC_AND_FILTERING.out.versions.ifEmpty([]) )

    // // Mixing bulk and single cell channels for clonal analysis
    // ch_repertoires_for_clones = ch_bulk_filtered
    //                                 .mix(SINGLE_CELL_QC_AND_FILTERING.out.repertoires)
    //                                 .dump(tag: 'sc bulk mix')

    // // Clonal analysis
    // CLONAL_ANALYSIS(
    //     ch_repertoires_for_clones,
    //     VDJ_ANNOTATION.out.imgt.collect(),
    //     ch_report_logo_img.collect().ifEmpty([])
    // )
    // ch_versions = ch_versions.mix( CLONAL_ANALYSIS.out.versions.ifEmpty([]))

    // if (!params.skip_report){
    //     REPERTOIRE_ANALYSIS_REPORTING(
    //         ch_presto_filterseq_logs.collect().ifEmpty([]),
    //         ch_presto_maskprimers_logs.collect().ifEmpty([]),
    //         ch_presto_pairseq_logs.collect().ifEmpty([]),
    //         ch_presto_clustersets_logs.collect().ifEmpty([]),
    //         ch_presto_buildconsensus_logs.collect().ifEmpty([]),
    //         ch_presto_postconsensus_pairseq_logs.collect().ifEmpty([]),
    //         ch_presto_assemblepairs_logs.collect().ifEmpty([]),
    //         ch_presto_collapseseq_logs.collect().ifEmpty([]),
    //         ch_presto_splitseq_logs.collect().ifEmpty([]),
    //         ch_reassign_logs.collect().ifEmpty([]),
    //         VDJ_ANNOTATION.out.changeo_makedb_logs.collect().ifEmpty([]),
    //         VDJ_ANNOTATION.out.logs.collect().ifEmpty([]),
    //         BULK_QC_AND_FILTER.out.logs.collect().ifEmpty([]),
    //         SINGLE_CELL_QC_AND_FILTERING.out.logs.collect().ifEmpty([]),
    //         CLONAL_ANALYSIS.out.logs.collect().ifEmpty([]),
    //         CLONAL_ANALYSIS.out.repertoire,
    //         ch_input.collect(),
    //         ch_report_rmd.collect(),
    //         ch_report_css.collect(),
    //         ch_report_logo.collect(),
    //         ch_validated_samplesheet.collect()
    //     )
    // }
    // ch_versions = ch_versions.mix( REPERTOIRE_ANALYSIS_REPORTING.out.versions )
    // ch_versions.dump(tag: "channel_versions")
    // // Software versions
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    // if (!params.skip_multiqc) {
    //     workflow_summary    = WorkflowAirrflow.paramsSummaryMultiqc(workflow, summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     methods_description    = WorkflowAirrflow.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    //     ch_methods_description = Channel.value(methods_description)

    //     ch_multiqc_files = Channel.empty()
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //     ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_html.ifEmpty([]))
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_fastp_json.ifEmpty([]))
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_postassembly_mqc.collect{it[1]}.ifEmpty([]))

    //     MULTIQC (
    //         ch_multiqc_files.collect(),
    //         ch_multiqc_config.collect(),
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         ch_report_logo.collect().ifEmpty([])
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)

    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
