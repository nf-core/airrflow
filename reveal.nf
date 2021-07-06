#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/bcellmagic
========================================================================================
 nf-core/bcellmagic Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/bcellmagic
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
checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Please provide input file containing the sample metadata with the '--input' option." }

// Input validation
if (params.input)  { 
      file(params.input, checkIfExists: true) 
} else { 
      exit 1, "Missing mandatory input: --input." 
}
 
if (params.miairr)  { 
      file(params.miairr, checkIfExists: true) 
} 

// If paths to DBS are provided 
if (params.igblast_base){
    Channel.fromPath("${params.igblast_base}")
    .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
    .set { ch_igblast }
}
if (params.imgtdb_base){
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

// Local: Sub-workflows
include { REVEAL_INPUT_CHECK } from './subworkflows/reveal_input_check'       addParams( options: [:] )

// Modules: local
include { GET_SOFTWARE_VERSIONS     } from './modules/local/get_software_versions'  addParams( options: [publish_files : ['csv':'']] )
include { IMMCANTATION  } from './modules/local/reveal/immcantation_container_version' addParams( options: [:] )
include { CHANGEO_CONVERTDB_FASTA } from './modules/local/changeo/changeo_convertdb_fasta'  addParams( options: modules['changeo_convertdb_fasta_from_airr'] )
include { FETCH_DATABASES } from './modules/local/fetch_databases'              addParams( options: [:] )
include { CHANGEO_ASSIGNGENES_REVEAL } from './modules/local/reveal/changeo_assigngenes_reveal'      addParams( options: modules['changeo_assigngenes_reveal'] )
include { CHANGEO_MAKEDB } from './modules/local/changeo/changeo_makedb'                addParams( options: modules['changeo_makedb_reveal'] ) 
include { FILTER_QUALITY  } from './modules/local/reveal/filter_quality' addParams( options: modules['filter_quality_reveal'] )
include { CHANGEO_PARSEDB_SPLIT } from './modules/local/changeo/changeo_parsedb_split'  addParams( options: modules['changeo_parsedb_split_reveal'] )
include { FILTER_JUNCTION_MOD3  } from './modules/local/reveal/filter_junction_mod3' addParams( options: modules['filter_quality_reveal'] )
include { CHIMERIC  } from './modules/local/reveal/chimeric' addParams( options: modules['filter_quality_reveal'] )
include { ADD_META_TO_TAB  } from './modules/local/reveal/add_meta_to_tab' addParams( options: modules['filter_quality_reveal'] )
include { COLLAPSE_DUPLICATES  } from './modules/local/reveal/collapse_duplicates' addParams( options: modules['filter_quality_reveal'] )



// include { CHANGEO_ASSIGNGENES } from './modules/local/changeo/changeo_assign_genes'  addParams( options: modules['changeo_assign_genes'] )

// nf-core/modules: Modules
include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options )


workflow REVEAL {

    ch_software_versions = Channel.empty()
    
    if (params.custom_container) {
        IMMCANTATION()
        ch_software_versions = ch_software_versions.mix(IMMCANTATION.out.version.first().ifEmpty(null))
    }

    // SUBWORKFLOW: Read in samplesheet, validate
    // and emit channels for fasta and tsv files
    REVEAL_INPUT_CHECK (ch_input, params.miairr, params.collapseby, params.cloneby)

    // If reassign requested, generate fasta from the tsv files    
    if (params.reassign) {
        ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA(REVEAL_INPUT_CHECK.out.ch_tsv).fasta
        ch_software_versions = ch_software_versions.mix(CHANGEO_CONVERTDB_FASTA.out.version.first().ifEmpty(null))
    } else {
        ch_fasta_from_tsv = Channel.empty()
    }

    // mix all fasta
    ch_fasta = REVEAL_INPUT_CHECK.out.ch_fasta.mix(ch_fasta_from_tsv)

   // FETCH DATABASES
   // TODO: this can take a long time, and the progress shows 0%. Would be
   // nice to have some better progress reporting.
   // And maybe run this as 2 separate steps, one for IMGT and one for IgBLAST?
    if (!params.igblast_base | !params.imgtdb_base) {
        FETCH_DATABASES()
        ch_software_versions = ch_software_versions.mix(FETCH_DATABASES.out.version.first().ifEmpty(null))
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
    }
    
    // Run Igblast for gene assignment
    CHANGEO_ASSIGNGENES_REVEAL (
        ch_fasta,
        ch_igblast.collect()
    )
    //ch_software_versions = ch_software_versions.mix(CHANGEO_ASSIGNGENES_REVEAL.out.version.first().ifEmpty(null))
    
    // Parse IgBlast results
    CHANGEO_MAKEDB (
        CHANGEO_ASSIGNGENES_REVEAL.out.fasta,
        CHANGEO_ASSIGNGENES_REVEAL.out.blast,
        ch_imgt.collect()
    )

    // Apply quality filters
    FILTER_QUALITY(CHANGEO_MAKEDB.out.tab)

    // Select only productive sequences and
    // sequences with junction length multiple of 3
    if (params.productive_only) {
        CHANGEO_PARSEDB_SPLIT (
            FILTER_QUALITY.out.tab
        )
        ch_repertoire = FILTER_JUNCTION_MOD3(CHANGEO_PARSEDB_SPLIT.out.tab).tab
    } else {
        ch_repertoire = FILTER_QUALITY.out.tab
    }

    // For bulk datasets, remove chimeric sequences
    // if requested
    if (params.remove_chimeric) {
        ch_repertoire
        .branch { it ->
            single: it[0].single_cell == 'true'
            bulk:   it[0].single_cell == 'false'
        }
        .set{ch_repertoire_by_processing}
        CHIMERIC(
            ch_repertoire_by_processing.bulk,
            ch_imgt.collect()
        ).tab
        // Mix with single
        ch_chimeric_pass = ch_repertoire_by_processing.single.mix(CHIMERIC.out.tab)
    } else {
        ch_chimeric_pass = ch_repertoire
    }

    ch_annotated_repertoires = ADD_META_TO_TAB(ch_chimeric_pass, REVEAL_INPUT_CHECK.out.validated_input)

    // Collapse duplicates
    // https://www.nextflow.io/docs/latest/operator.html#grouptuple
    ch_collapsable = ch_annotated_repertoires.tab
     .map{ it -> [ it[0].single_cell, it[0], it[1] ] }
     .groupTuple(by: [0])
     .map{ it -> [it[1], it[2].toList()] }
     .dump()

    COLLAPSE_DUPLICATES(ch_collapsable,params.collapseby)    
    

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
        //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        
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