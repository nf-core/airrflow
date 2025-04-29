include { PREPARE_TRUST4_REFERENCE                                      } from '../../modules/local/prepare_trust4_reference'
include { TRUST4                                                        } from '../../modules/nf-core/trust4/main'
include { FASTQ_INPUT_CHECK                                             } from '../../subworkflows/local/fastq_input_check'
include { CHANGEO_PARSEDB_SELECT_LOCUS                                  } from '../../modules/local/changeo/changeo_parsedb_select_locus'
include { CHANGEO_CONVERTDB_FASTA as CHANGEO_CONVERTDB_FASTA_FROM_AIRR  } from '../../modules/local/changeo/changeo_convertdb_fasta'
include { FASTP                                                         } from '../../modules/nf-core/fastp/main'
include { RENAME_FASTQ as RENAME_FASTQ_TRUST4                           } from '../../modules/local/rename_fastq'



workflow RNASEQ_INPUT {

    take:
    ch_input
    ch_igblast_reference

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


    // validate library generation method parameters
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
        error "The TRUST4 library generation method does not require this reference, please provide a compliant reference file instead or select another library method option."
    }

    // Fastp
    save_merged = false
    FASTP (
        ch_reads,
        [],
        [],
        save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    ch_rename_fastq = FASTP.out.reads.map { meta, reads -> [meta, reads[0], reads[1]] }

    // rename fastp output
    RENAME_FASTQ_TRUST4(
        ch_rename_fastq
    )

    ch_reads_fastp_filtered = RENAME_FASTQ_TRUST4.out.reads.dump(tag: "fastp_filtered")

    PREPARE_TRUST4_REFERENCE(
        ch_reads_fastp_filtered.first(),
        ch_igblast_reference
    )


    // create trust4 input
    ch_reads_trust4 = ch_reads_fastp_filtered.map{ meta, read_1, read_2  -> [ meta, [], [read_1, read_2] ] }

    PREPARE_TRUST4_REFERENCE.out.trust4_reference.dump(tag: "trust4_reference")

    ch_reads_trust4.dump(tag: "trust4_input")

    // create barcode and umi channels for nf-core trust4 module
    ch_reference_trust4 = PREPARE_TRUST4_REFERENCE.out.trust4_reference.map { reference -> [[id: "igblast_reference"], reference] }
    barcode_channel = ch_reads_fastp_filtered.map { meta, read_1, read_2 ->  [meta, params.trust4_cell_barcode_read] }
    umi_channel = ch_reads_fastp_filtered.map { meta, read_1, read_2 -> [meta, params.trust4_umi_read] }
    barcodewhitelist_channel = ch_reads_fastp_filtered.map { meta, read_1, read_2 -> 
                                            [meta, params.trust4_barcode_whitelist ? file(params.trust4_barcode_whitelist) : []] }   

    TRUST4(
        ch_reads_trust4,
        ch_reference_trust4.collect(),
        Channel.of([[], []]).collect(),
        barcode_channel,
        umi_channel,
        barcodewhitelist_channel
    )

    ch_trust4_out = TRUST4.out.outs

    // check whether input is sc or bulk and extract respective airr file for downstream processing
    ch_trust4_out
        .branch {
            meta, out_files ->
                bulk : meta["single_cell"] == "false"
                    return [ meta, out_files.find { it.endsWith("${meta.id}_airr.tsv") } ]
                sc : meta["single_cell"] == "true"
                    return [ meta, out_files.find { it.endsWith("${meta.id}_barcode_airr.tsv") } ]
        }
        .set { ch_trust4_airr_file }


    // create channel with airr file
    ch_trust4_airr_file.bulk.mix ( ch_trust4_airr_file.sc ).set { ch_trust4_airr }

    // select only provided locus
    CHANGEO_PARSEDB_SELECT_LOCUS(ch_trust4_airr)


    // convert airr tsv to fasta
    CHANGEO_CONVERTDB_FASTA_FROM_AIRR(
                CHANGEO_PARSEDB_SELECT_LOCUS.out.tab
            )

    ch_fasta = CHANGEO_CONVERTDB_FASTA_FROM_AIRR.out.fasta


    emit:
    versions = ch_versions
    // fastp
    fastp_reads_json = FASTP.out.json.collect{ meta,json -> json }
    fastp_reads_html = FASTP.out.html.collect{ meta,html -> html }
    // complete trust4 output
    outs = ch_trust4_out
    // trust4 airr file
    airr = ch_trust4_airr
    // trust4 output converted to FASTA format
    fasta = ch_fasta
    samplesheet = FASTQ_INPUT_CHECK.out.samplesheet

}
