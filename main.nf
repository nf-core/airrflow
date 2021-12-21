#!/usr/bin/env nextflow
/*
========================================================================================
                        nf-core/bcellmagic
========================================================================================
    GitHub  : https://github.com/nf-core/bcellmagic
    Website : https://nf-co.re/bcellmagic
    Slack   : https://nfcore.slack.com/channels/bcellmagic
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
include { REVEAL } from './workflows/reveal'

workflow {
    if (params.subworkflow == "bcellmagic") {
        include { BCELLMAGIC } from './workflows/bcellmagic'
        BCELLMAGIC()
    } else if (params.subworkflow == "reveal") {
        REVEAL()
    } else {
        exit 1
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
