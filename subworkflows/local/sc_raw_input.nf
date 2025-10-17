include { CELLRANGER_VDJ                                                } from '../../modules/nf-core/cellranger/vdj/main'
include { UNZIP_CELLRANGERDB                                            } from '../../modules/local/unzip_cellrangerdb'
include { RENAME_FILE as RENAME_FILE_TSV                                } from '../../modules/local/rename_file'
include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR  } from '../../modules/local/changeo/changeo_convertdb_fasta'
include { FASTQ_INPUT_CHECK                                             } from '../../subworkflows/local/fastq_input_check'


workflow SC_RAW_INPUT {

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

    // validate library generation method parameter
    if (params.vprimers) {
        error "The single-cell 10X genomics library generation method does not require V-region primers, please provide a reference file instead or select another library method option."
    } else if (params.race_linker) {
        error "The single-cell 10X genomics library generation method does not require the --race_linker parameter, please provide a reference file instead or select another library method option."
    }
    if (params.cprimers)  {
        error "The single-cell 10X genomics library generation method does not require C-region primers, please provide a reference file instead or select another library method option."
    }
    if (params.umi_length > 0)  {
        error "The single-cell 10X genomics library generation method does not require to set the UMI length, please provide a reference file instead or select another library method option."
    }
    if (params.reference_10x)  {
        // necessary to allow tar.gz files as input so that tests can run
        if (params.reference_10x.endsWith(".tar.gz")){
            UNZIP_CELLRANGERDB(
                params.reference_10x
            )
            ch_versions = ch_versions.mix(UNZIP_CELLRANGERDB.out.versions)
            UNZIP_CELLRANGERDB.out.unzipped.set { ch_sc_reference }
        } else {
            ch_sc_reference = Channel.fromPath(params.reference_10x, checkIfExists: true)
        }
    } else {
        error "The single-cell 10X genomics library generation method requires you to provide a reference file."
    }

    // run cellranger vdj
    CELLRANGER_VDJ (
        ch_reads,
        ch_sc_reference.collect()
    )
    ch_versions = ch_versions.mix(CELLRANGER_VDJ.out.versions)

    ch_cellranger_out = CELLRANGER_VDJ.out.outs

    ch_cellranger_out
        .map { meta, out_files ->
                [ meta, out_files.find { it.endsWith("airr_rearrangement.tsv") } ]
            }
        .set { ch_cellranger_airr }

    // TODO : add VALIDATE_INPUT Module
    // this module requires input in csv format... Might need to create this in an extra module

    // rename tsv file to unique name
    RENAME_FILE_TSV(
                ch_cellranger_airr
            )
        .set { ch_renamed_tsv }

    // convert airr tsv to fasta (cellranger does not create any fasta with clonotype information)
    CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                RENAME_FILE_TSV.out.file
            )

    ch_versions = ch_versions.mix(CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.versions)

    ch_fasta = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta

    // TODO: here you can add support for MiXCR sc protocols.


    emit:
    versions = ch_versions
    // complete cellranger output
    outs = ch_cellranger_out
    // cellranger output in airr format
    airr = ch_cellranger_airr
    // cellranger output converted to FASTA format
    fasta = ch_fasta
    samplesheet = FASTQ_INPUT_CHECK.out.samplesheet
}
