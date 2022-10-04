#!/usr/bin/env nextflow
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/


if (params.miairr)  {
    file(params.miairr, checkIfExists: true)
}

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config  = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_report_logo = Channel.fromPath(params.report_logo, checkIfExists: true)


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Modules: local
include { IMMCANTATION  } from '../modules/local/reveal/immcantation_container_version'
include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR } from '../modules/local/changeo/changeo_convertdb_fasta'
include { FETCH_DATABASES } from '../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_IMGT } from '../modules/local/unzip_db'
include { CHANGEO_ASSIGNGENES_REVEAL } from '../modules/local/reveal/changeo_assigngenes_reveal'
include { CHANGEO_MAKEDB_REVEAL } from '../modules/local/reveal/changeo_makedb_reveal'
include { FILTER_QUALITY  } from '../modules/local/reveal/filter_quality'
include { CHANGEO_PARSEDB_SPLIT as CHANGEO_PARSEDB_SPLIT_REVEAL} from '../modules/local/changeo/changeo_parsedb_split'
include { FILTER_JUNCTION_MOD3  } from '../modules/local/reveal/filter_junction_mod3'
include { CHANGEO_CREATEGERMLINES_REVEAL as CREATEGERMLINES } from '../modules/local/reveal/changeo_creategermlines_reveal'
include { REMOVE_CHIMERIC  } from '../modules/local/enchantr/remove_chimeric'
include { SINGLE_CELL_QC  } from '../modules/local/enchantr/single_cell_qc'
include { ADD_META_TO_TAB  } from '../modules/local/reveal/add_meta_to_tab'
include { DETECT_CONTAMINATION  } from '../modules/local/enchantr/detect_contamination'
include { COLLAPSE_DUPLICATES  } from '../modules/local/enchantr/collapse_duplicates'
include { FIND_THRESHOLD } from '../modules/local/enchantr/find_threshold'
include { DEFINE_CLONES } from '../modules/local/enchantr/define_clones'
include { DOWSER_LINEAGES } from '../modules/local/enchantr/dowser_lineages'
include { REPORT_FILE_SIZE     } from '../modules/local/enchantr/report_file_size'

// Local: Sub-workflows
include { REVEAL_INPUT_CHECK } from '../subworkflows/local/reveal_input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow REVEAL {

    log.warn "\n----------\nREVEAL lifecycle stage: experimental.\n----------\n"

    ch_versions = Channel.empty()
    ch_file_sizes = Channel.empty()

    IMMCANTATION()
    ch_versions = ch_versions.mix(IMMCANTATION.out.versions)

    // SUBWORKFLOW: Read in samplesheet, validate
    // and emit channels for fasta and tsv files
    REVEAL_INPUT_CHECK (ch_input, params.miairr, params.collapseby, params.cloneby, params.reassign)

    // If reassign requested, generate fasta from the tsv files
    if (params.reassign) {
        CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
            REVEAL_INPUT_CHECK.out.ch_tsv
        )
        ch_fasta_from_tsv = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta
        ch_versions = ch_versions.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.versions.ifEmpty(null))
        ch_file_sizes = ch_file_sizes.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.logs)
    } else {
        ch_fasta_from_tsv = Channel.empty()
    }


    // mix all fasta
    ch_fasta = REVEAL_INPUT_CHECK.out.ch_fasta.mix(ch_fasta_from_tsv)

    // FETCH DATABASES
    // TODO: this can take a long time, and the progress shows 0%. Would be
    // nice to have some better progress reporting.
    // And maybe run this as 2 separate steps, one for IMGT and one for IgBLAST?
    if( params.igblast_base ){

        if (params.igblast_base.endsWith(".zip")) {
            Channel.fromPath("${params.igblast_base}")
                    .ifEmpty{ exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions.ifEmpty(null))
        } else {
            Channel.fromPath("${params.igblast_base}")
                .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
                .set { ch_igblast }
        }
    }
    if( params.imgtdb_base ){

        if (params.imgtdb_base.endsWith(".zip")) {
            Channel.fromPath("${params.imgtdb_base}")
                    .ifEmpty{ exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
                    .set { ch_imgt_zipped }
            UNZIP_IMGT( ch_imgt_zipped.collect() )
            ch_imgt = UNZIP_IMGT.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IMGT.out.versions.ifEmpty(null))
        } else {
            Channel.fromPath("${params.imgtdb_base}")
                .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
                .set { ch_imgt }
        }
    }

    if (!params.igblast_base | !params.imgtdb_base) {
        FETCH_DATABASES()
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_imgt = FETCH_DATABASES.out.imgt
        ch_versions = ch_versions.mix(FETCH_DATABASES.out.versions.ifEmpty(null))
    }

    // Run Igblast for gene assignment
    CHANGEO_ASSIGNGENES_REVEAL (
        ch_fasta,
        ch_igblast.collect()
    )
    ch_file_sizes = ch_file_sizes.mix(CHANGEO_ASSIGNGENES_REVEAL.out.logs)
    ch_versions = ch_versions.mix(CHANGEO_ASSIGNGENES_REVEAL.out.versions.ifEmpty(null))

    // Parse IgBlast results
    CHANGEO_MAKEDB_REVEAL (
        CHANGEO_ASSIGNGENES_REVEAL.out.fasta,
        CHANGEO_ASSIGNGENES_REVEAL.out.blast,
        ch_imgt.collect()
    )
    ch_file_sizes = ch_file_sizes.mix(CHANGEO_MAKEDB_REVEAL.out.logs)
    ch_versions = ch_versions.mix(CHANGEO_MAKEDB_REVEAL.out.versions.ifEmpty(null))

    // Apply quality filters
    // TODO: mv to enchantr and emit versions
    FILTER_QUALITY(CHANGEO_MAKEDB_REVEAL.out.tab)
    ch_file_sizes = ch_file_sizes.mix(FILTER_QUALITY.out.logs)

    // Select only productive sequences and
    // sequences with junction length multiple of 3
    if (params.productive_only) {
        CHANGEO_PARSEDB_SPLIT_REVEAL (
            FILTER_QUALITY.out.tab
        )
        ch_file_sizes = ch_file_sizes.mix(CHANGEO_PARSEDB_SPLIT_REVEAL.out.logs)
        ch_versions = ch_versions.mix(CHANGEO_PARSEDB_SPLIT_REVEAL.out.versions.ifEmpty(null))

        // TODO: Add to enchantr and emit versions?
        FILTER_JUNCTION_MOD3(
            CHANGEO_PARSEDB_SPLIT_REVEAL.out.tab
        )
        ch_file_sizes = ch_file_sizes.mix(FILTER_JUNCTION_MOD3.out.logs)
        ch_repertoire = FILTER_JUNCTION_MOD3.out.tab.ifEmpty(null)
    } else {
        ch_repertoire = FILTER_QUALITY.out.tab.ifEmpty(null)
    }

    // Add metadata to the rearrangement files, to be used later
    // for grouping, subsetting, plotting....
    ADD_META_TO_TAB(
        ch_repertoire,
        REVEAL_INPUT_CHECK.out.validated_input.collect()
    )
    ch_file_sizes = ch_file_sizes.mix(ADD_META_TO_TAB.out.logs)

    ch_repertoire_by_processing = ADD_META_TO_TAB.out.tab
        .dump(tag: 'meta_to_tab_out')
        .branch { it ->
            single: it[0].single_cell == 'true'
            bulk:   it[0].single_cell == 'false'
        }

    ch_repertoire_by_processing.bulk
    .dump(tag: 'bulk')

    ch_repertoire_by_processing.single
    .dump(tag: 'single')

    // For bulk datasets, remove chimeric sequences
    // if requested
    if (params.remove_chimeric) {

        // Create germlines (not --cloned)
        CREATEGERMLINES(
            ch_repertoire_by_processing.bulk,
            ch_imgt.collect()
        )
        ch_file_sizes = ch_file_sizes.mix(CREATEGERMLINES.out.logs)
        ch_versions = ch_versions.mix(CREATEGERMLINES.out.versions.ifEmpty(null))

        // Remove chimera
        REMOVE_CHIMERIC(
            CREATEGERMLINES.out.tab,
            ch_imgt.collect()
        )
        ch_file_sizes = ch_file_sizes.mix(REMOVE_CHIMERIC.out.logs)
        ch_bulk_chimeric_pass = REMOVE_CHIMERIC.out.tab
        ch_versions = ch_versions.mix(REMOVE_CHIMERIC.out.versions.ifEmpty(null))

    } else {
        ch_bulk_chimeric_pass = ch_repertoire_by_processing.bulk
    }

    // For Bulk data, detect cross-contamination
    // This is only informative at this time
    // TODO: add a flag to specify remove suspicious sequences
    // and update file size log accordingly
    DETECT_CONTAMINATION(
        ch_bulk_chimeric_pass
            .map{ it -> [ it[1] ] }
            .collect(),
        'id')
    // TODO file size
    ch_versions = ch_versions.mix(DETECT_CONTAMINATION.out.versions.ifEmpty(null))

    COLLAPSE_DUPLICATES(
        ch_bulk_chimeric_pass
            .map{ it -> [ it[1] ] }
            .collect(),
        params.collapseby
    )
    ch_versions = ch_versions.mix(COLLAPSE_DUPLICATES.out.versions.ifEmpty(null))
    // TODO file size
    // TODO channel by params.cloneby


    // For single cell, specific QC
    // analyze all files together, looking for overlaps
    SINGLE_CELL_QC(
        ch_repertoire_by_processing.single
        .map{ it -> [ it[1] ] }
        .collect()
    )
    ch_file_sizes = ch_file_sizes.mix(SINGLE_CELL_QC.out.logs)
    ch_versions = ch_versions.mix(SINGLE_CELL_QC.out.versions.ifEmpty(null))

    // If params.threshold is auto,
    // 1) use distToNearest and findThreshold to determine
    // the threshold that will be used to identify sets of clonally
    // related sequences. If threshold found, continue, to 2), else,
    // stop and report a threshold could be identified.
    // 2) create a report with plots of the distToNearest distribution
    // and the threshold.
    // Else
    // Use the threshold to find clones, grouping by params.cloneby and
    // create a report

    if (params.threshold == "auto") {
        FIND_THRESHOLD (
            COLLAPSE_DUPLICATES.out.tab.mix(SINGLE_CELL_QC.out.tab)
            .map{ it -> [ it[1] ] }
            .collect(),
            params.cloneby,
            params.singlecell
        )
        ch_threshold = FIND_THRESHOLD.out.mean_threshold

        clone_threshold = ch_threshold
            .splitText( limit:1 ) { it.trim().toString() }
            .dump(tag: 'clone_threshold')
            .filter { it != 'NA'}
            .dump(tag: "threshold")
            .ifEmpty { exit 1, "Automatic clone_threshold is 'NA'. Consider setting params.threshold manually."}

    } else {
        // TODO: Fix * --threshold: expected type: String, found: Integer (1)
        clone_threshold = params.threshold
    }

    DEFINE_CLONES(
        SINGLE_CELL_QC.out.tab.mix(COLLAPSE_DUPLICATES.out.tab)
        .map{ it -> [ it[1] ] }
        .collect(),
        params.cloneby,
        params.singlecell,
        clone_threshold,
        ch_imgt.collect()
    )

    DOWSER_LINEAGES(
        DEFINE_CLONES.out.tab
        .map{ it -> [ it[1] ] }
        .flatten()
    )

    // TODO fix file sizes
    // Process logs to report file sizes at each step
    //REPORT_FILE_SIZE (
    //    ch_file_sizes.map { it }.collect()
    //)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowBcellmagic.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //   ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.yaml.collect())

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.collect(),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_report_logo.collect().ifEmpty([])
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
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
