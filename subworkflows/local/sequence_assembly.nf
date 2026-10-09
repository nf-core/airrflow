#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Local: Sub-workflows
include { FASTQ_INPUT_CHECK           } from '../../subworkflows/local/fastq_input_check'
include { PRESTO_UMI                  } from '../../subworkflows/local/presto_umi'
include { PRESTO_SANS_UMI             } from '../../subworkflows/local/presto_sans_umi'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQUENCE_ASSEMBLY {

    take:
    ch_input // channel: reads
    ch_igblast
    library_generation_method
    adapter_fasta
    maskprimers_extract
    vprimers
    cprimers
    race_linker
    umi_length
    internal_cregion_sequences
    maskprimers_align_race
    index_file
    umi_position
    umi_start
    collapseby
    cloneby
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
    filterseq_q
    buildconsensus_maxerror
    buildconsensus_maxgap
    primer_consensus

    main:

    // Validate params
    if (!library_generation_method) {
        error "Please specify a library generation method with the `--library_generation_method` option."
    }

    if (adapter_fasta) {
        ch_adapter_fasta = channel.fromPath(adapter_fasta, checkIfExists: true)
    } else {
        ch_adapter_fasta = []
    }


    // Validate library generation method parameter
    if (library_generation_method == 'specific_pcr_umi'){
        if (!maskprimers_extract){
            if (vprimers)  {
                ch_vprimers_fasta = channel.fromPath(vprimers, checkIfExists: true)
            } else {
                error "Please provide a V-region primers fasta file with the '--vprimers' option when using the 'specific_pcr_umi' library generation method."
            }
            if (cprimers)  {
                ch_cprimers_fasta = channel.fromPath(cprimers, checkIfExists: true)
            } else {
                error "Please provide a C-region primers fasta file with the '--cprimers' option when using the 'specific_pcr_umi' library generation method."
            }
            if (race_linker)  {
                error "Please do not set '--race_linker' when using the 'specific_pcr_umi' library generation method."
            }
        } else {
            ch_vprimers_fasta = channel.of([])
            ch_cprimers_fasta = channel.of([])
        }
        if (umi_length < 2)  {
            error "The 'specific_pcr_umi' library generation method requires setting the '--umi_length' to a value greater than 1."
        }
        if (internal_cregion_sequences) {
            ch_internal_cregion = channel.fromPath(internal_cregion_sequences, checkIfExists: true)
        } else {
            ch_internal_cregion = channel.of([])
        }
    } else if (library_generation_method == 'specific_pcr') {
        if (vprimers)  {
            ch_vprimers_fasta = channel.fromPath(vprimers, checkIfExists: true)
        } else {
            error "Please provide a V-region primers fasta file with the '--vprimers' option when using the 'specific_pcr' library generation method."
        }
        if (cprimers)  {
            ch_cprimers_fasta = channel.fromPath(cprimers, checkIfExists: true)
        } else {
            error "Please provide a C-region primers fasta file with the '--cprimers' option when using the 'specific_pcr' library generation method."
        }
        if (race_linker)  {
            error "Please do not set '--race_linker' when using the 'specific_pcr' library generation method."
        }
        if (umi_length > 0)  {
            error "Please do not set a UMI length with the library preparation method 'specific_pcr'. Please specify instead a method that suports umi."
        } else {
            umi_length = 0
        }
        if (internal_cregion_sequences) {
            error "Please do not set '--internal_cregion_sequences' when using the 'specific_pcr' library generation method without UMIs."
        }
    } else if (library_generation_method == 'dt_5p_race_umi') {
        if (vprimers) {
            error "The oligo-dT 5'-RACE UMI library generation method does not accept V-region primers, please provide a linker with '--race_linker' instead or select another library method option."
        } else if (race_linker) {
            ch_vprimers_fasta = channel.fromPath(race_linker, checkIfExists: true)
        } else if (maskprimers_align_race) {
            ch_vprimers_fasta = channel.of([])
        } else {
            error "The oligo-dT 5'-RACE UMI library generation method requires a linker or Template Switch Oligo sequence, please provide it with the option '--race_linker'."
        }
        if (cprimers)  {
            ch_cprimers_fasta = channel.fromPath(cprimers, checkIfExists: true)
        } else {
            error "The oligo-dT 5'-RACE UMI library generation method requires the C-region primer sequences, please provide a fasta file with the '--cprimers' option."
        }
        if (umi_length < 2)  {
            error "The oligo-dT 5'-RACE UMI 'dt_5p_race_umi' library generation method requires specifying the '--umi_length' to a value greater than 1."
        }
        if (internal_cregion_sequences) {
            ch_internal_cregion = channel.fromPath(internal_cregion_sequences, checkIfExists: true)
        } else {
            ch_internal_cregion = channel.of([])
        }
    } else if (library_generation_method == 'dt_5p_race') {
        if (vprimers) {
            error "The oligo-dT 5'-RACE library generation method does not accept V-region primers, please provide a linker with '--race_linker' instead or select another library method option."
        } else if (race_linker) {
            ch_vprimers_fasta = channel.fromPath(race_linker, checkIfExists: true)
        } else if (maskprimers_align_race) {
            ch_vprimers_fasta = channel.of([])
        } else {
            error "The oligo-dT 5'-RACE library generation method requires a linker or Template Switch Oligo sequence, please provide it with the option '--race_linker'."
        }
        if (cprimers)  {
            ch_cprimers_fasta = channel.fromPath(cprimers, checkIfExists: true)
        } else {
            error "The oligo-dT 5'-RACE library generation method requires the C-region primer sequences, please provide a fasta file with the '--cprimers' option."
        }
        if (umi_length > 0)  {
            error "Please do not set a UMI length with the library preparation method oligo-dT 5'-RACE 'dt_5p_race'. Please specify instead a method that suports umi (e.g. 'dt_5p_race_umi')."
        } else {
            umi_length = 0
        }
        if (internal_cregion_sequences) {
            error "Please do not set '--internal_cregion_sequences' when using the 'dt_5p_race' library generation method without UMIs."
        }
    } else {
        error "The provided library generation method is not supported. Please check the docs for `--library_generation_method`."
    }

    // Validate UMI position
    if (index_file & umi_position == 'R2') {error "Please do not set `--umi_position` option if index file with UMIs is provided."}
    if (umi_length < 0) {error "Please provide the UMI barcode length in the option `--umi_length`. To run without UMIs, set umi_length to 0."}
    if (!index_file & umi_start != 0) {error "Setting a UMI start position is only allowed when providing the UMIs in a separate index read file. If so, please provide the `--index_file` flag as well."}

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_versions = channel.empty()

    FASTQ_INPUT_CHECK(
        ch_input,
        library_generation_method,
        collapseby,
        cloneby,
        index_file
    )
    ch_versions = ch_versions.mix(FASTQ_INPUT_CHECK.out.versions)

    ch_reads = FASTQ_INPUT_CHECK.out.reads

    if (umi_length == 0) {
        //
        // SUBWORKFLOW: pRESTO without UMIs
        //
        PRESTO_SANS_UMI(
            ch_reads,
            ch_cprimers_fasta,
            ch_vprimers_fasta,
            ch_adapter_fasta,
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
            filterseq_q
        )
        ch_presto_fasta = PRESTO_SANS_UMI.out.fasta
        ch_fastp_reads_html = PRESTO_SANS_UMI.out.fastp_reads_html
        ch_fastp_reads_json = PRESTO_SANS_UMI.out.fastp_reads_json
        ch_fastqc_postassembly = PRESTO_SANS_UMI.out.fastqc_postassembly_gz
        ch_presto_assemblepairs_logs = PRESTO_SANS_UMI.out.presto_assemblepairs_logs
        ch_presto_filterseq_logs = PRESTO_SANS_UMI.out.presto_filterseq_logs
        ch_presto_maskprimers_logs = PRESTO_SANS_UMI.out.presto_maskprimers_logs
        ch_presto_collapseseq_logs = PRESTO_SANS_UMI.out.presto_collapseseq_logs
        ch_presto_splitseq_logs = PRESTO_SANS_UMI.out.presto_splitseq_logs
        ch_presto_pairseq_logs = channel.empty()
        ch_presto_clustersets_logs = channel.empty()
        ch_presto_buildconsensus_logs = channel.empty()
        ch_presto_postconsensus_pairseq_logs = channel.empty()
        ch_presto_UMIreads = channel.empty()

    } else {
        //
        // SUBWORKFLOW: pRESTO with UMIs
        //
        PRESTO_UMI (
            ch_reads,
            ch_cprimers_fasta,
            ch_vprimers_fasta,
            ch_adapter_fasta,
            ch_internal_cregion,
            ch_igblast.collect(),
            maskprimers_align_race,
            umi_position,
            cprimer_position,
            index_file,
            save_trimmed,
            primer_maxlen,
            primer_r1_maxerror,
            primer_r1_mask_mode,
            umi_length,
            umi_start,
            primer_r2_extract_len,
            primer_r2_mask_mode,
            maskprimers_align,
            primer_r2_maxerror,
            primer_revpr,
            maskprimers_extract,
            cprimer_start,
            primer_r1_extract_len,
            vprimer_start,
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
        ch_presto_fasta = PRESTO_UMI.out.fasta
        ch_fastp_reads_html = PRESTO_UMI.out.fastp_reads_html
        ch_fastp_reads_json = PRESTO_UMI.out.fastp_reads_json
        ch_fastqc_postassembly = PRESTO_UMI.out.fastqc_postassembly_gz
        ch_presto_filterseq_logs = PRESTO_UMI.out.presto_filterseq_logs
        ch_presto_maskprimers_logs = PRESTO_UMI.out.presto_maskprimers_logs
        ch_presto_pairseq_logs = PRESTO_UMI.out.presto_pairseq_logs
        ch_presto_clustersets_logs = PRESTO_UMI.out.presto_clustersets_logs
        ch_presto_buildconsensus_logs = PRESTO_UMI.out.presto_buildconsensus_logs
        ch_presto_postconsensus_pairseq_logs = PRESTO_UMI.out.presto_postconsensus_pairseq_logs
        ch_presto_assemblepairs_logs = PRESTO_UMI.out.presto_assemblepairs_logs
        ch_presto_collapseseq_logs = PRESTO_UMI.out.presto_collapseseq_logs
        ch_presto_splitseq_logs = PRESTO_UMI.out.presto_splitseq_logs
        ch_presto_UMIreads = PRESTO_UMI.out.presto_UMIreads
    }


    emit:
    // assembled sequences in fasta format
    fasta = ch_presto_fasta
    // validated metadata
    samplesheet = FASTQ_INPUT_CHECK.out.samplesheet
    //fastp
    fastp_reads_html = ch_fastp_reads_html
    fastp_reads_json = ch_fastp_reads_json
    // fastqc files for multiQC report
    fastqc_postassembly = ch_fastqc_postassembly
    // presto logs for html report
    presto_filterseq_logs = ch_presto_filterseq_logs
    presto_maskprimers_logs = ch_presto_maskprimers_logs
    presto_pairseq_logs = ch_presto_pairseq_logs
    presto_clustersets_logs = ch_presto_clustersets_logs
    presto_buildconsensus_logs = ch_presto_buildconsensus_logs
    presto_postconsensus_pairseq_logs = ch_presto_postconsensus_pairseq_logs
    presto_assemblepairs_logs = ch_presto_assemblepairs_logs
    presto_collapseseq_logs = ch_presto_collapseseq_logs
    presto_splitseq_logs = ch_presto_splitseq_logs
    presto_UMIreads = ch_presto_UMIreads
    versions = ch_versions
}
