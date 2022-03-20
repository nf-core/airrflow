#!/usr/bin/env nextflow
/*
========================================================================================
                        nf-core/airrflow
========================================================================================
    GitHub  : https://github.com/nf-core/airrflow
    Website : https://nf-co.re/airrflow
    Slack   : https://nfcore.slack.com/channels/airrflow
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

if (params.subworkflow == 'bcellmagic') {
    include { BCELLMAGIC } from './workflows/bcellmagic'
} else if (params.subworkflow == 'reveal') {
    include { REVEAL } from './workflows/reveal'
}

workflow NFCORE_AIRRFLOW {
    if (params.subworkflow == "bcellmagic") {
        BCELLMAGIC()
    } else if (params.subworkflow == "reveal") {
        REVEAL()
    } else {
        exit 1
    }
}

workflow {
    NFCORE_AIRRFLOW()
}

/*
========================================================================================
    THE END
========================================================================================
*/
