include { TRUST4                                                        } from '../../modules/local/trust4'
include { FASTQ_INPUT_CHECK                                             } from '../../subworkflows/local/fastq_input_check'
include { RENAME_FILE as RENAME_FILE_TSV                                } from '../../modules/local/rename_file'
include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR  } from '../../modules/local/changeo/changeo_convertdb_fasta'

workflow RNASEQ_INPUT {

    take:
    ch_input

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    //
    // read in samplesheet, validate and stage input fies
    //
    FASTQ_INPUT_CHECK(
        ch_input
    )
    ch_versions = ch_versions.mix(FASTQ_INPUT_CHECK.out.versions)

    ch_reads = FASTQ_INPUT_CHECK.out.reads

    // validate library generation method parameteç
    if (params.vprimers) {
        error "The TRUST4 library generation method does not require V-region primers, please provide a reference file instead or select another library method option."
    } else if (params.race_linker) {
        error "The TRUST4 10X genomics library generation method does not require the --race_linker parameter, please provide a reference file instead or select another library method option."
    }
    if (params.cprimers)  {
        error "The TRUST4 library generation method does not require C-region primers, please provide a reference file instead or select another library method option."
    }
    if (params.umi_length > 0)  {
        error "TRUST4 library generation method does not require to set the UMI length, please provide a reference file instead or select another library method option."
    }
    if (params.reference_10x)  {
        // necessary to allow tar.gz files as input so that tests can run
        error "The TRUST4 library generation method does not require this reference, please provide a compliant reference file instead or select another library method option."
    }
    if (!params.coord_fasta) {
        error "Please provide a reference file for the TRUST4 library generation method."
    }


    ch_reads.map{ meta, input_file  ->
        [ meta, [], input_file ] }
    .set { ch_reads_new }

    // ch_reads_new.view()

    // if (params.vdj_reference != null) {
    //     ch_vdjref = Channel.of([[], file(params.vdj_reference)])
    // }
    // else {
    //     ch_vdjref = Channel.of([[], []])
    // }

    TRUST4(
        ch_reads_new,
        Channel.of([[], file(params.coord_fasta)]),
        Channel.of([[], []])
    )

    ch_trust4_out = TRUST4.out.out
    // ch_trust4_airr = TRUST4.out.airr_tsv
    ch_trust4_airr = TRUST4.out.barcode_airr

    ch_trust4_airr.view()
    // rename tsv file to unique name
    RENAME_FILE_TSV(
                ch_trust4_airr
            )
        .set { ch_renamed_tsv }

    // convert airr tsv to fasta (cellranger does not create any fasta with clonotype information)
    CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                RENAME_FILE_TSV.out.file
            )

    ch_fasta = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta

    emit:
    versions = ch_versions
    // complete trust4 output
    outs = ch_trust4_out
    // trust4 airr file
    airr = ch_trust4_airr
    // trust4 output converted to FASTA format
    fasta = ch_fasta
    samplesheet = FASTQ_INPUT_CHECK.out.samplesheet
}
