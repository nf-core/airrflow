#!/usr/bin/env nextflow
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBcellmagic.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Please provide input file containing the sample metadata with the '--input' option." }

if (!params.library_generation_method) {
    exit 1, "Please specify a library generation method with the `--library_generation_method` option."
}

// Validate library generation method parameter
if (params.library_generation_method == 'specific_pcr_umi'){
    if (params.vprimers)  {
        ch_vprimers_fasta = Channel.fromPath(params.vprimers, checkIfExists: true)
    } else {
        exit 1, "Please provide a V-region primers fasta file with the '--vprimers' option when using the 'specific_pcr_umi' library generation method."
    }
    if (params.cprimers)  {
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true)
    } else {
        exit 1, "Please provide a C-region primers fasta file with the '--cprimers' option when using the 'specific_pcr_umi' library generation method."
    }
    if (params.race_linker)  {
        exit 1, "Please do not set '--race_linker' when using the 'specific_pcr_umi' library generation method."
    }
    if (params.umi_length < 2)  {
        exit 1, "The 'specific_pcr_umi' library generation method requires setting the '--umi_length' to a value greater than 1."
    }
} else if (params.library_generation_method == 'specific_pcr') {
    if (params.vprimers)  {
        ch_vprimers_fasta = Channel.fromPath(params.vprimers, checkIfExists: true)
    } else {
        exit 1, "Please provide a V-region primers fasta file with the '--vprimers' option when using the 'specific_pcr' library generation method."
    }
    if (params.cprimers)  {
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true)
    } else {
        exit 1, "Please provide a C-region primers fasta file with the '--cprimers' option when using the 'specific_pcr' library generation method."
    }
    if (params.race_linker)  {
        exit 1, "Please do not set '--race_linker' when using the 'specific_pcr' library generation method."
    }
    if (params.umi_length > 0)  {
        exit 1, "Please do not set a UMI length with the library preparation method 'specific_pcr'. Please specify instead a method that suports umi."
    } else {
        params.umi_length = 0
    }

} else if (params.library_generation_method == 'dt_5p_race_umi') {
    if (params.vprimers) {
        exit 1, "The oligo-dT 5'-RACE UMI library generation method does not accept V-region primers, please provide a linker with '--race_linker' instead or select another library method option."
    } else if (params.race_linker) {
        ch_vprimers_fasta = Channel.fromPath(params.race_linker, checkIfExists: true)
    } else {
        exit 1, "The oligo-dT 5'-RACE UMI library generation method requires a linker or Template Switch Oligo sequence, please provide it with the option '--race_linker'."
    }
    if (params.cprimers)  {
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true)
    } else {
        exit 1, "The oligo-dT 5'-RACE UMI library generation method requires the C-region primer sequences, please provide a fasta file with the '--cprimers' option."
    }
    if (params.umi_length < 2)  {
        exit 1, "The oligo-dT 5'-RACE UMI 'dt_5p_race_umi' library generation method requires specifying the '--umi_length' to a value greater than 1."
    }
} else if (params.library_generation_method == 'dt_5p_race') {
    if (params.vprimers) {
        exit 1, "The oligo-dT 5'-RACE library generation method does not accept V-region primers, please provide a linker with '--race_linker' instead or select another library method option."
    } else if (params.race_linker) {
        ch_vprimers_fasta = Channel.fromPath(params.race_linker, checkIfExists: true)
    } else {
        exit 1, "The oligo-dT 5'-RACE library generation method requires a linker or Template Switch Oligo sequence, please provide it with the option '--race_linker'."
    }
    if (params.cprimers)  {
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true)
    } else {
        exit 1, "The oligo-dT 5'-RACE library generation method requires the C-region primer sequences, please provide a fasta file with the '--cprimers' option."
    }
    if (params.umi_length > 0)  {
        exit 1, "Please do not set a UMI length with the library preparation method oligo-dT 5'-RACE 'dt_5p_race'. Please specify instead a method that suports umi (e.g. 'dt_5p_race_umi')."
    } else {
        params.umi_length = 0
    }
} else {
    exit 1, "The provided library generation method is not supported. Please check the docs for `--library_generation_method`."
}

// Validate UMI position
if (params.index_file & params.umi_position == 'R2') {exit 1, "Please do not set `--umi_position` option if index file with UMIs is provided."}
if (params.umi_length < 0) {exit 1, "Please provide the UMI barcode length in the option `--umi_length`. To run without UMIs, set umi_length to 0."}
if (!params.index_file & params.umi_start != 0) {exit 1, "Setting a UMI start position is only allowed when providing the UMIs in a separate index read file. If so, please provide the `--index_file` flag as well."}

// If paths to databases are provided
if( params.igblast_base ){
    Channel.fromPath("${params.igblast_base}")
            .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
            .set { ch_igblast }
}
if( params.imgtdb_base ){
    Channel.fromPath("${params.imgtdb_base}")
            .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
            .set { ch_imgt }
}


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Rmarkdown report file
ch_rmarkdown_report      = Channel.fromPath( ["$projectDir/assets/repertoire_comparison.Rmd",
                                    "$projectDir/assets/references.bibtex",
                                    "$projectDir/assets/nf-core_style.css",
                                    "$projectDir/assets/nf-core-bcellmagic_logo.png"],
                                    checkIfExists: true).dump(tag: 'report files')

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'  addParams( options: [publish_files : ['csv':'']] )

//CHANGEO
include { FETCH_DATABASES } from '../modules/local/fetch_databases'              addParams( options: [:] )
include { CHANGEO_ASSIGNGENES } from '../modules/local/changeo/changeo_assigngenes'      addParams( options: modules['changeo_assigngenes'] )
include { CHANGEO_MAKEDB } from '../modules/local/changeo/changeo_makedb'                addParams( options: modules['changeo_makedb'] )
include { CHANGEO_PARSEDB_SPLIT } from '../modules/local/changeo/changeo_parsedb_split'  addParams( options: modules['changeo_parsedb_split'] )
include { CHANGEO_PARSEDB_SELECT } from '../modules/local/changeo/changeo_parsedb_select'    addParams( options: modules['changeo_parsedb_select'] )
include { CHANGEO_CONVERTDB_FASTA } from '../modules/local/changeo/changeo_convertdb_fasta'  addParams( options: modules['changeo_convertdb_fasta'] )

//SHAZAM
include { SHAZAM_TIGGER_THRESHOLD } from '../modules/local/shazam/shazam_tigger_threshold'  addParams( options: modules['shazam_tigger_threshold'] )

//CHANGEO
include { CHANGEO_DEFINECLONES } from '../modules/local/changeo/changeo_defineclones'        addParams( options: modules['changeo_defineclones'] )
include { CHANGEO_CREATEGERMLINES } from '../modules/local/changeo/changeo_creategermlines'  addParams( options: modules['changeo_creategermlines'] )
include { CHANGEO_BUILDTREES } from '../modules/local/changeo/changeo_buildtrees'        addParams( options: modules['changeo_buildtrees'] )

//ALAKAZAM
include { ALAKAZAM_LINEAGE } from '../modules/local/alakazam/alakazam_lineage'            addParams( options: modules['alakazam_lineage'] )
include { ALAKAZAM_SHAZAM_REPERTOIRES } from '../modules/local/alakazam/alakazam_shazam_repertoires'   addParams ( options: modules['alakazam_shazam_repertoires'] )

//LOG PARSING
include { PARSE_LOGS } from '../modules/local/parse_logs'                           addParams( options: modules['parse_logs'] )

// Local: Sub-workflows
include { INPUT_CHECK           } from '../subworkflows/local/input_check'          addParams( options: [:] )
include { MERGE_TABLES_WF       } from '../subworkflows/local/merge_tables_wf'      addParams( options: modules['merge_tables'] )
include { PRESTO_UMI            } from '../subworkflows/local/presto_umi'           addParams( options: [:] )
include { PRESTO_SANS_UMI            } from '../subworkflows/local/presto_sans_umi'           addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                } from '../modules/nf-core/modules/fastqc/main'        addParams( options: modules['fastqc'] )
include { MULTIQC               } from '../modules/nf-core/modules/multiqc/main'       addParams( options: multiqc_options )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BCELLMAGIC {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( ch_input )
        .groupTuple(by: [0])
        .map{ it -> [ it[0], it[1].flatten() ] }
        .set{ ch_fastqc }

    ch_presto = ch_fastqc.map{ it -> it.flatten() }

    //
    // MODULE: FastQC
    //
    FASTQC ( ch_fastqc )

    // Channel for software versions
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    if (params.umi_length == 0) {
        //
        // SUBWORKFLOW: pRESTO without UMIs
        //
        PRESTO_SANS_UMI (
            ch_presto,
            ch_cprimers_fasta,
            ch_vprimers_fasta
        )
        ch_presto_fasta = PRESTO_SANS_UMI.out.fasta
        ch_presto_software = PRESTO_SANS_UMI.out.software
        ch_fastqc_postassembly_gz = PRESTO_SANS_UMI.out.fastqc_postassembly_gz
        ch_presto_assemblepairs_logs = PRESTO_SANS_UMI.out.presto_assemblepairs_logs
        ch_presto_filterseq_logs = PRESTO_SANS_UMI.out.presto_filterseq_logs
        ch_presto_maskprimers_logs = PRESTO_SANS_UMI.out.presto_maskprimers_logs
        ch_presto_collapseseq_logs = PRESTO_SANS_UMI.out.presto_collapseseq_logs
        ch_presto_splitseq_logs = PRESTO_SANS_UMI.out.presto_splitseq_logs
        // These channels will be empty in the sans-UMI workflow
        ch_presto_pairseq_logs = Channel.empty()
        ch_presto_clustersets_logs = Channel.empty()
        ch_presto_buildconsensus_logs = Channel.empty()
        ch_presto_postconsensus_pairseq_logs = Channel.empty()
    } else {
        //
        // SUBWORKFLOW: pRESTO with UMIs
        //
        PRESTO_UMI (
            ch_presto,
            ch_cprimers_fasta,
            ch_vprimers_fasta
        )
        ch_presto_fasta = PRESTO_UMI.out.fasta
        ch_presto_software = PRESTO_UMI.out.software
        ch_fastqc_postassembly_gz = PRESTO_UMI.out.fastqc_postassembly_gz
        ch_presto_filterseq_logs = PRESTO_UMI.out.presto_filterseq_logs
        ch_presto_maskprimers_logs = PRESTO_UMI.out.presto_maskprimers_logs
        ch_presto_pairseq_logs = PRESTO_UMI.out.presto_pairseq_logs
        ch_presto_clustersets_logs = PRESTO_UMI.out.presto_clustersets_logs
        ch_presto_buildconsensus_logs = PRESTO_UMI.out.presto_buildconsensus_logs
        ch_presto_postconsensus_pairseq_logs = PRESTO_UMI.out.presto_postconsensus_pairseq_logs
        ch_presto_assemblepairs_logs = PRESTO_UMI.out.presto_assemblepairs_logs
        ch_presto_collapseseq_logs = PRESTO_UMI.out.presto_collapseseq_logs
        ch_presto_splitseq_logs = PRESTO_UMI.out.presto_splitseq_logs
    }

    ch_software_versions = ch_software_versions.mix(ch_presto_software)

    // FETCH DATABASES
    if (!params.igblast_base | !params.imgtdb_base) {
        FETCH_DATABASES()
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
        ch_software_versions = ch_software_versions.mix(FETCH_DATABASES.out.version.first().ifEmpty(null))
    }

    // Run Igblast for gene assignment
    CHANGEO_ASSIGNGENES (
        ch_presto_fasta,
        ch_igblast.collect()
    )
    ch_software_versions = ch_software_versions.mix(CHANGEO_ASSIGNGENES.out.version.first().ifEmpty(null))

    // Make IgBlast results table
    CHANGEO_MAKEDB (
        CHANGEO_ASSIGNGENES.out.fasta,
        CHANGEO_ASSIGNGENES.out.blast,
        ch_imgt.collect()
    )

    // Select only productive sequences.
    CHANGEO_PARSEDB_SPLIT (
        CHANGEO_MAKEDB.out.tab
    )

    // Selecting IGH for ig loci, TR for tr loci.
    CHANGEO_PARSEDB_SELECT(
        CHANGEO_PARSEDB_SPLIT.out.tab
    )

    // Convert sequence table to fasta.
    CHANGEO_CONVERTDB_FASTA (
        CHANGEO_PARSEDB_SELECT.out.tab
    )

    // Subworkflow: merge tables from the same patient
    MERGE_TABLES_WF(CHANGEO_PARSEDB_SELECT.out.tab)

    // Shazam clonal threshold and tigger genotyping
    SHAZAM_TIGGER_THRESHOLD(
        MERGE_TABLES_WF.out.tab.dump(tag: 'merge tables output'),
        ch_imgt.collect()
    )

    ch_software_versions = ch_software_versions.mix(SHAZAM_TIGGER_THRESHOLD.out.version.first().ifEmpty(null)).dump()

    // Define B-cell clones
    CHANGEO_DEFINECLONES(
        SHAZAM_TIGGER_THRESHOLD.out.tab,
        SHAZAM_TIGGER_THRESHOLD.out.threshold,
        SHAZAM_TIGGER_THRESHOLD.out.fasta
    )

    // Identify germline sequences
    CHANGEO_CREATEGERMLINES(
        CHANGEO_DEFINECLONES.out.tab,
        CHANGEO_DEFINECLONES.out.fasta,
        ch_imgt.collect()
    )

    // Lineage reconstruction alakazam
    if (!params.skip_lineage) {
        ALAKAZAM_LINEAGE(
            CHANGEO_CREATEGERMLINES.out.tab.dump(tag:'creategermlines_output')
        )
    }

    ch_software_versions = ch_software_versions.mix(ALAKAZAM_LINEAGE.out.version.first().ifEmpty(null)).dump()

    ch_all_tabs_repertoire = CHANGEO_CREATEGERMLINES.out.tab
                                                    .map{ it -> [ it[1] ] }
                                                    .collect()
                                                    .dump(tag:'repertoire_all')

    // Process logs parsing: getting sequence numbers
    PARSE_LOGS(
        ch_presto_filterseq_logs.collect(),
        ch_presto_maskprimers_logs.collect(),
        ch_presto_pairseq_logs.collect().ifEmpty([]),
        ch_presto_clustersets_logs.collect().ifEmpty([]),
        ch_presto_buildconsensus_logs.collect().ifEmpty([]),
        ch_presto_postconsensus_pairseq_logs.collect().ifEmpty([]),
        ch_presto_assemblepairs_logs.collect(),
        ch_presto_collapseseq_logs.collect(),
        ch_presto_splitseq_logs.collect(),
        CHANGEO_MAKEDB.out.logs.collect(),
        CHANGEO_DEFINECLONES.out.logs.collect(),
        CHANGEO_CREATEGERMLINES.out.logs.collect(),
        ch_input
    )

    // Alakazam shazam repertoire comparison report
    if (!params.skip_report){
        ALAKAZAM_SHAZAM_REPERTOIRES(
            ch_all_tabs_repertoire,
            PARSE_LOGS.out.logs.collect(),
            ch_rmarkdown_report.collect()
        )
    }

    // Software versions
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowBcellmagic.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_config)
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_postassembly_gz.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
