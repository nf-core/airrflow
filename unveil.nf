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

igblast_db = params.igblast_db
imgt_db = params.imgt_db
igblastn = params.igblastn

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
include { UNVEIL_INPUT_CHECK } from './subworkflows/unveil_input_check'       addParams( options: [:] )

// Modules: local
include { GET_SOFTWARE_VERSIONS     } from './modules/local/get_software_versions'  addParams( options: [publish_files : ['csv':'']] )
include { IMMCANTATION  } from './modules/local/unveil/immcantation_container_version' addParams( options: [:] )
include { CHANGEO_CONVERTDB_FASTA } from './modules/local/changeo/changeo_convertdb_fasta'  addParams( options: modules['changeo_convertdb_fasta_from_airr'] )
// include { CHANGEO_ASSIGNGENES } from './modules/local/changeo/changeo_assign_genes'  addParams( options: modules['changeo_assign_genes'] )

// nf-core/modules: Modules
include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options )


workflow UNVEIL {

    ch_software_versions = Channel.empty()
    
    if (workflow.container) {
        IMMCANTATION()
        ch_software_versions = ch_software_versions.mix(IMMCANTATION.out.version.first().ifEmpty(null))
    }

    // SUBWORKFLOW: Read in samplesheet, validate
    // and emit channels for fasta and tsv files
    UNVEIL_INPUT_CHECK (ch_input, params.miairr, params.collapseby, params.cloneby)

    // If reassign requested, generate fasta from the tsv files    
    if (params.reassign) {
        ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA(UNVEIL_INPUT_CHECK.out.ch_tsv)
        ch_software_versions = ch_software_versions.mix(CHANGEO_CONVERTDB_FASTA.out.version.first().ifEmpty(null))
    } else {
        ch_fasta_from_tsv = Channel.empty()
    }

    // mix all fasta
    ch_fasta = UNVEIL_INPUT_CHECK.out.ch_fasta.mix(ch_fasta_from_tsv)
   
    // Assign genes
    //CHANGEO_ASSIGNGENES (ch_fasta, igblast_db, imgt_db, igblastn)
    //ch_software_versions = ch_software_versions.mix(CHANGEO_ASSIGNGENES.out.version.first().ifEmpty(null))
    
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