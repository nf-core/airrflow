#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/bcellmagic
========================================================================================
Documentation: https://nf-co.re/bcellmagic
Code: https://github.com/nf-core/bcellmagic
----------------------------------------------------------------------------------------
*/

////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Please provide input file containing the sample metadata with the '--input' option." }

// Validate primer protocol
if (params.protocol == "pcr_umi"){
    if (params.vprimers)  { 
        ch_vprimers_fasta = Channel.fromPath(params.vprimers, checkIfExists: true) 
    } else { 
        exit 1, "Please provide a V-region primers fasta file with the '--vprimers' option, or specify a 5'RACE protocol with the '--protocol' option." 
    }
    if (params.cprimers)  { 
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true) 
    } else { 
        exit 1, "Please provide a C-region primers fasta file with the '--cprimers' option." 
    }
} else if (params.protocol == "race_5p_umi") {
    if (params.vprimers) { 
        exit 1, "The 5' RACE protocol does not accept V-region primers, please remove the option '--vprimers' or provide another protocol."
    } else if (params.race_linker) {
        ch_vprimers_fasta = Channel.fromPath(params.race_linker, checkIfExists: true)
    } else {
        exit 1, "The 5' RACE protocol requires a linker or Template Switch Oligo sequence, please provide it with the option '--race_linker'."
    }
    if (params.cprimers)  { 
        ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true) 
    } else { 
        exit 1, "Please provide a C-region primers fasta file with the '--cprimers' option." 
    }
}

// Validate UMI position
if (params.index_file & params.umi_position == 'R2') {exit 1, "Please do not set `--umi_position` option if index file with UMIs is provided."}
if (params.umi_length == 0) {exit 1, "Please provide the UMI barcode length in the option `--umi_length`."}
if (!params.index_file & params.umi_start != 0) {exit 1, "Setting a UMI start position is only allowed when providing the UMIs in a separate index read file. If so, please provide the `--index_file` flag as well."}

// If paths to DBS are provided 
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


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions'  addParams( options: [publish_files : ['csv':'']] )
include { MERGE_UMI } from './modules/local/merge_UMI'                          addParams( options: [:] )
include { RENAME_FASTQ } from './modules/local/rename_fastq'                    addParams( options: [:] )
include { GUNZIP } from './modules/local/gunzip'                                addParams( options: [:] )
//PRESTO
include { PRESTO_FILTERSEQ } from './modules/local/presto/presto_filterseq'            addParams( options: modules['presto_filterseq'] )
include { PRESTO_MASKPRIMERS } from './modules/local/presto/presto_maskprimers'        addParams( options: modules['presto_maskprimers'] )
include { PRESTO_PAIRSEQ } from './modules/local/presto/presto_pairseq'                addParams( options: modules['presto_pairseq'] )
include { PRESTO_CLUSTERSETS } from './modules/local/presto/presto_clustersets'        addParams( options: modules['presto_clustersets'] )
include { PRESTO_PARSE_CLUSTER } from './modules/local/presto/presto_parse_cluster'    addParams( options: [:] )
include { PRESTO_BUILDCONSENSUS } from './modules/local/presto/presto_buildconsensus'  addParams( options: modules['presto_buildconsensus'] )
include { PRESTO_POSTCONSENSUS_PAIRSEQ } from './modules/local/presto/presto_postconsensus_pairseq'    addParams( options: modules['presto_postconsensus_pairseq'] )
include { PRESTO_ASSEMBLEPAIRS } from './modules/local/presto/presto_assemblepairs'    addParams( options: modules['presto_assemblepairs'] )
include { PRESTO_PARSEHEADERS as PRESTO_PARSEHEADERS_COLLAPSE } from './modules/local/presto/presto_parseheaders'  addParams( options: modules['presto_parseheaders_collapse'] )
include { PRESTO_PARSEHEADERS_PRIMERS } from './modules/local/presto/presto_parseheaders_primers'      addParams( options: [:] )
include { PRESTO_PARSEHEADERS_METADATA } from './modules/local/presto/presto_parseheaders_metadata'    addParams( options: [:] )
include { PRESTO_COLLAPSESEQ } from './modules/local/presto/presto_collapseseq'        addParams( options: modules['presto_collapseseq'] )
include { PRESTO_SPLITSEQ } from './modules/local/presto/presto_splitseq'              addParams( options: modules['presto_splitseq'] )
//CHANGEO
include { FETCH_DATABASES } from './modules/local/fetch_databases'              addParams( options: [:] )
include { CHANGEO_ASSIGNGENES } from './modules/local/changeo/changeo_assigngenes'      addParams( options: modules['changeo_assigngenes'] )
include { CHANGEO_MAKEDB } from './modules/local/changeo/changeo_makedb'                addParams( options: modules['changeo_makedb'] ) 
include { CHANGEO_PARSEDB_SPLIT } from './modules/local/changeo/changeo_parsedb_split'  addParams( options: modules['changeo_parsedb_split'] )
include { CHANGEO_PARSEDB_SELECT } from './modules/local/changeo/changeo_parsedb_select'    addParams( options: modules['changeo_parsedb_select'] )
include { CHANGEO_CONVERTDB_FASTA } from './modules/local/changeo/changeo_convertdb_fasta'  addParams( options: modules['changeo_convertdb_fasta'] )
//SHAZAM
include { SHAZAM_TIGGER_THRESHOLD } from './modules/local/shazam/shazam_tigger_threshold'  addParams( options: modules['shazam_tigger_threshold'] )
//CHANGEO
include { CHANGEO_DEFINECLONES } from './modules/local/changeo/changeo_defineclones'        addParams( options: modules['changeo_defineclones'] )
include { CHANGEO_CREATEGERMLINES } from './modules/local/changeo/changeo_creategermlines'  addParams( options: modules['changeo_creategermlines'] )
include { CHANGEO_BUILDTREES } from './modules/local/changeo/changeo_buildtrees'        addParams( options: modules['changeo_buildtrees'] )
//ALAKAZAM
include { ALAKAZAM_LINEAGE } from './modules/local/alakazam/alakazam_lineage'            addParams( options: modules['alakazam_lineage'] )
include { ALAKAZAM_SHAZAM_REPERTOIRES } from './modules/local/alakazam/alakazam_shazam_repertoires'   addParams ( options: modules['alakazam_shazam_repertoires'] )
//LOG PARSING
include { PARSE_LOGS } from './modules/local/parse_logs'                        addParams( options: modules['parse_logs'] )

// Local: Sub-workflows
include { INPUT_CHECK           } from './subworkflows/input_check'       addParams( options: [:] )
include { MERGE_TABLES_WF       } from './subworkflows/merge_tables_wf'      addParams( options: modules['merge_tables'] )

// nf-core/modules: Modules
include { FASTQC                } from './modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc'] )
include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options )

////////////////////////////////////////////////////
/*           BCELLMAGIC WORKFLOW                  */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow BCELLMAGIC {

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input )
    .groupTuple(by: [0])
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_fastqc }

    ch_merge_umi_gunzip = ch_fastqc.map{ it -> it.flatten() }

    // FastQC
    FASTQC ( ch_fastqc )

    // Channel for software versions
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    
    // Merge UMI from index file to R1 if provided
    if (params.index_file) {
        MERGE_UMI ( ch_merge_umi_gunzip )
        .set{ ch_gunzip }
    } else {
        RENAME_FASTQ ( ch_merge_umi_gunzip )
        .set{ ch_gunzip }
    }

    // gunzip fastq.gz to fastq
    GUNZIP ( ch_gunzip )
    ch_software_versions = ch_software_versions.mix(GUNZIP.out.version.first().ifEmpty(null))

    // Filter sequences by quality score
    PRESTO_FILTERSEQ ( GUNZIP.out.reads )
    ch_software_versions = ch_software_versions.mix(PRESTO_FILTERSEQ.out.version.first().ifEmpty(null))

    // Mask primers
    PRESTO_MASKPRIMERS ( 
        PRESTO_FILTERSEQ.out.reads,
        ch_cprimers_fasta.collect(),
        ch_vprimers_fasta.collect()
    )

    // Pre-consensus pair
    PRESTO_PAIRSEQ (
        PRESTO_MASKPRIMERS.out.reads
    )

    // Cluster sequences by similarity
    PRESTO_CLUSTERSETS (
        PRESTO_PAIRSEQ.out.reads
    )
    // ch_software_versions = ch_software_versions.mix(PRESTO_CLUSTERSETS.out.version.first().ifEmpty(null))

    // Annotate cluster into barcode field
    PRESTO_PARSE_CLUSTER (
        PRESTO_CLUSTERSETS.out.reads
    )

    // Build consensus of sequences with same UMI barcode
    PRESTO_BUILDCONSENSUS (
        PRESTO_PARSE_CLUSTER.out.reads
    )

    // Post-consensus pair 
    PRESTO_POSTCONSENSUS_PAIRSEQ (
        PRESTO_BUILDCONSENSUS.out.reads
    )

    // Assemble read pairs
    PRESTO_ASSEMBLEPAIRS (
        PRESTO_POSTCONSENSUS_PAIRSEQ.out.reads
    )

    // Combine UMI duplicate count
    PRESTO_PARSEHEADERS_COLLAPSE (
        PRESTO_ASSEMBLEPAIRS.out.reads
    )

    // Annotate primers in C_PRIMER and V_PRIMER field
    PRESTO_PARSEHEADERS_PRIMERS (
        PRESTO_PARSEHEADERS_COLLAPSE.out.reads
    )

    // Annotate metadata on primer headers
    PRESTO_PARSEHEADERS_METADATA (
        PRESTO_PARSEHEADERS_PRIMERS.out.reads
    )

    // Mark and count duplicate sequences with different UMI barcodes (DUPCOUNT)
    PRESTO_COLLAPSESEQ (
        PRESTO_PARSEHEADERS_METADATA.out.reads
    )

    // Filter out sequences with less than 2 representative duplicates with different UMIs
    PRESTO_SPLITSEQ (
        PRESTO_COLLAPSESEQ.out.reads
    )

    // FETCH DATABASES
    if (!params.igblast_base | !params.imgtdb_base) {
        FETCH_DATABASES()
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
    }

    // Run Igblast for gene assignment
    CHANGEO_ASSIGNGENES (
        PRESTO_SPLITSEQ.out.fasta,
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
    CHANGEO_PARSEDB_SELECT (
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
        MERGE_TABLES_WF.out,
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

    //Changeo build trees
    //CHANGEO_BUILDTREES(
    //    CHANGEO_CREATEGERMLINES.out.tab
    //)

    // Lineage reconstruction alakazam
    ALAKAZAM_LINEAGE(
        CHANGEO_CREATEGERMLINES.out.tab
    )

    ch_software_versions = ch_software_versions.mix(ALAKAZAM_LINEAGE.out.version.first().ifEmpty(null)).dump()

    //ch_all_tabs_repertoire = ALAKAZAM_LINEAGE.out.tab.collect()

    // Alakazam shazam repertoire comparison
    //ALAKAZAM_SHAZAM_REPERTOIRES(
    //    ch_all_tabs_repertoire
    //)

    // Process logs parsing: getting sequence numbers
    PARSE_LOGS(
        PRESTO_FILTERSEQ.out.logs.collect(),
        PRESTO_MASKPRIMERS.out.logs.collect(),
        PRESTO_PAIRSEQ.out.logs.collect(),
        PRESTO_CLUSTERSETS.out.logs.collect(),
        PRESTO_BUILDCONSENSUS.out.logs.collect(),
        PRESTO_POSTCONSENSUS_PAIRSEQ.out.logs.collect(),
        PRESTO_ASSEMBLEPAIRS.out.logs.collect(),
        PRESTO_COLLAPSESEQ.out.logs.collect(),
        PRESTO_SPLITSEQ.out.logs.collect(),
        CHANGEO_MAKEDB.out.logs.collect(),
        CHANGEO_DEFINECLONES.out.logs.collect(),
        CHANGEO_CREATEGERMLINES.out.logs.collect(),
        ch_input
    )

    // Software versions
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        
        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////