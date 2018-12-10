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


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/bcellmagic v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/bcellmagic --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --cprimers                    Path to CPrimers FASTA File
      --vprimers                    Path to VPrimers FASTA File
      --metadata                    Path to Metadata TSV
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --imgtdb_base                 Path to predownloaded imgtDB 
      --igblast_base                Path to predownloaded igblastDB

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
output_docs = Channel.fromPath("$baseDir/docs/output.md")
//Defaults for igblast
igblast_base = false
imgtdb_base = false

// If paths to DBS are provided 
if( params.igblast_base ){
    Channel.fromPath("${params.igblast_base}")
    .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
    .set { ch_igblast_db_for_process_igblast_mix }
}
if( params.imgtdb_base ){
    Channel.fromPath("${params.imgtdb_base}")
    .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
    .into { ch_imgt_db_for_igblast_filter_mix;ch_imgt_db_for_shazam_mix;ch_imgt_db_for_germline_sequences_mix }
}
//Set up channels for input primers
Channel.fromPath(file("${params.cprimers}"))
       .ifEmpty{exit 1, "Please specify CPRimers FastA File!"}
       .set {ch_cprimers_fasta}
Channel.fromPath(file("${params.vprimers}"))
       .ifEmpty{exit 1, "Please specify VPrimers FastA File!"}
       .set { ch_vprimers_fasta }

saveDBs = false

//Other parameters
filterseq_q=20


//Download data process
process fetchDBs{
    tag "fetchBlastDBs"

    publishDir path: { params.saveDBs ? "${params.outdir}/dbs" : params.outdir },
    saveAs: { params.saveDBs ? it : null }, mode: 'copy'

    when:
    !params.igblast_base | !params.imgtdb_base

    output:
    file "igblast_base" into ch_igblast_db_for_process_igblast
    file "imgtdb_base" into (ch_imgt_db_for_igblast_filter,ch_imgt_db_for_shazam,ch_imgt_db_for_germline_sequences)
    
    script:
    """
    fetch_igblastdb.sh -o igblast_base
    fetch_imgtdb.sh -o imgtdb_base
    imgt2igblast.sh -i imgtdb_base -o igblast_base
    """
}

/*
 * Create a channel for metadata and raw files
 * Columns = id, source, treatment, extraction_time, population, R1, R2, I1
 */
 file_meta = file(params.metadata)
 ch_read_files_for_merge_r1_umi = Channel.from(file_meta)
                .splitCsv(header: true, sep:'\t')
                .map { col -> tuple("${col.ID}", "${col.Source}", "${col.Treatment}","${col.Extraction_time}","${col.Population}",returnFile("${col.R1}"),returnFile("${col.R2}"),returnFile("${col.I1}"))}
                .dump()
// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/bcellmagic v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/bcellmagic'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['MetaData']        = params.metadata
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['IGDB Path']    = params.igblast_base
summary['IMGT Path']    = params.imgtdb_base
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

//Merge I1 UMIs into R1 file
process merge_r1_umi {
    tag "${R1.baseName}"

    input:
    set val(id), val(source), val(treatment), val(extraction_time), val(population), file(R1), file(R2), file(I1) from ch_read_files_for_merge_r1_umi

    output:
    file "*UMI_R1.fastq" into ch_fastqs_for_processing_umi
    file "${R2.baseName}" into ch_fastqs_for_processing_r2
    set val("$treatment"),val("$extraction_time"),val("$population") into ch_meta_env_for_anno
    val("$id") into (ch_sample_for_alakazam,ch_sample_for_shazam)

    script:
    """
    merge_R1_umi.py -R1 "${R1}" -I1 "${I1}" -o "${R1.baseName}_UMI_R1.fastq.gz"
    gunzip "${R1.baseName}_UMI_R1.fastq.gz"
    gunzip -f "${R2}"
    """
}


//Filter by Sequence Quality
process filter_by_sequence_quality {
    tag "${umi.baseName}"

    input:
    file(umi) from ch_fastqs_for_processing_umi
    file(r2) from ch_fastqs_for_processing_r2

    output:
    file "${umi.baseName}_quality-pass.fastq" into ch_filtered_by_seq_quality_for_primer_Masking_UMI
    file "${r2.baseName}_quality-pass.fastq" into ch_filtered_by_seq_quality_for_primerMasking_R2

    script:
    """
    FilterSeq.py quality -s $umi -q $filterseq_q --outname "${umi.baseName}"
    FilterSeq.py quality -s $r2 -q $filterseq_q --outname "${r2.baseName}"
    """
}

//Mask them primers
process mask_primers {
    tag "${umi_file.baseName}"

    input:
    file(umi_file) from ch_filtered_by_seq_quality_for_primer_Masking_UMI
    file(r2_file) from ch_filtered_by_seq_quality_for_primerMasking_R2
    file(cprimers) from ch_cprimers_fasta 
    file(vprimers) from ch_vprimers_fasta

    output:
    file "${umi_file.baseName}_UMI_R1_primers-pass.fastq" into ch_for_pair_seq_umi_file
    file "${r2_file.baseName}_R2_primers-pass.fastq" into ch_for_pair_seq_r2_file

    script:
    """
    MaskPrimers.py score -s $umi_file -p ${cprimers} --start 8 --mode cut --barcode --outname ${umi_file.baseName}_UMI_R1
    MaskPrimers.py score -s $r2_file -p ${vprimers} --start 0 --mode mask --outname ${r2_file.baseName}_R2
    """
}

//Pair the UMI_R1 and R2
process pair_seq{
    tag "${umi.baseName}"

    input:
    file umi from ch_for_pair_seq_umi_file
    file r2 from ch_for_pair_seq_r2_file

    output:
    file "${umi.baseName}_pair-pass.fastq" into ch_umi_for_umi_cluster_sets
    file "${r2.baseName}_pair-pass.fastq" into ch_r2_for_umi_cluster_sets

    script:
    """
    PairSeq.py -1 $umi -2 $r2 --1f BARCODE --coord illumina
    """
}

//Deal with too low UMI diversity
process cluster_sets {
    tag "${umi.baseName}"

    input:
    file umi from ch_umi_for_umi_cluster_sets
    file r2 from ch_r2_for_umi_cluster_sets

    output:
    file "${umi.baseName}_UMI_R1_cluster-pass.fastq" into ch_umi_for_reheader
    file "${r2.baseName}_R2_cluster-pass.fastq" into ch_r2_for_reheader

    script:
    """
    ClusterSets.py set -s $umi --outname ${umi.baseName}_UMI_R1 
    ClusterSets.py set -s $r2 --outname ${r2.baseName}_R2
    """
}

//ParseHeaders to annotate barcode into cluster names
process reheader {
    tag "${umi.baseName}" 

    input:
    file umi from ch_umi_for_reheader
    file r2 from ch_r2_for_reheader

    output:
    file "${umi.baseName}_reheader.fastq" into ch_umi_for_consensus
    file "${r2.baseName}_reheader.fastq" into ch_r2_for_consensus

    script:
    """
    ParseHeaders.py copy -s $umi -f BARCODE -k CLUSTER --act cat
    ParseHeaders.py copy -s $r2 -f BARCODE -k CLUSTER --act cat 
    """
}


//Build UMI consensus
process build_consensus{
    tag "${umi.baseName}"

    input:
    file umi from ch_umi_for_consensus
    file r2 from ch_r2_for_consensus

    output:
    file "${umi.baseName}_UMI_R1_consensus-pass.fastq" into ch_consensus_passed_umi
    file "${r2.baseName}_R2_consensus-pass.fastq" into ch_consensus_passed_r2

    script:
    """
    BuildConsensus.py -s $umi --bf CLUSTER --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${umi.baseName}_UMI_R1
    BuildConsensus.py -s $r2 --bf CLUSTER --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${r2.baseName}_R2
    """
}

//Repair again UMI_R1+R2
process repair{
    tag "${umi.baseName}"

    input:
    file umi from ch_consensus_passed_umi
    file r2 from ch_consensus_passed_r2

    output:
    file "*UMI_R1_consensus-pass_pair-pass.fastq" into ch_repaired_UMI_for_assembly
    file "*R2_consensus-pass_pair-pass.fastq" into ch_repaired_r2_for_assembly

    script:
    """
    PairSeq.py -1 $umi -2 $r2 --coord presto
    """
}

//Assemble the UMI consensus mate pairs
process assemble{
    tag "${umi.baseName}"

    input:
    file umi from ch_repaired_UMI_for_assembly
    file r2 from ch_repaired_r2_for_assembly

    output:
    file "${umi.baseName}_UMI_R1_R2_assemble-pass.fastq" into ch_for_combine_UMI

    script:
    """
    AssemblePairs.py align -1 $umi -2 $r2 --coord presto --rc tail --1f CONSCOUNT PRCONS --2f CONSCOUNT PRCONS --outname ${umi.baseName}_UMI_R1_R2
    """
}

//combine UMI read group size annotations
process combine_umi_read_groups{
    tag "${assembled.baseName}"

    input:
    file assembled from ch_for_combine_UMI

    output:
    file "*_reheader.fastq" into ch_for_prcons_parseheaders

    script:
    """
    ParseHeaders.py collapse -s $assembled -f CONSCOUNT --act min
    """
}

//Copy field PRCONS to have annotation for C_primer and V_primer independently
process copy_prcons{
    tag "${combined.baseName}"

    input:
    file combined from ch_for_prcons_parseheaders

    output:
    file "*reheader_reheader.fastq" into ch_for_metadata_anno

    script:
    """
    ParseHeaders.py copy -s $combined -f PRCONS PRCONS --act first last -k C_PRIMER V_PRIMER
    """
}

//Add Metadata annotation to headers
process metadata_anno{
    tag "${prcons.baseName}"

    input:
    file prcons from ch_for_metadata_anno
    set val(treatment),val(extraction_time),val(population) from ch_meta_env_for_anno

    output:
    file "*_reheader_reheader_reheader.fastq" into ch_for_dedup 

    script:
    """
    ParseHeaders.py add -s $prcons -f TREATMENT EXTRACT_TIME POPULATION -u $treatment $extraction_time $population
    """
}

//Removal of duplicate sequences
process dedup {
    tag "${dedup.baseName}"

    input:
    file dedup from ch_for_dedup

    output:
    file "${dedup.baseName}_UMI_R1_R2_collapse-unique.fastq" into ch_for_filtering

    script:
    """
    CollapseSeq.py -s $dedup -n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname ${dedup.baseName}_UMI_R1_R2
    """
}

//Filtering to sequences with at least two representative reads and convert to FastA
process filter_seqs{
    tag "${dedupped.baseName}"

    input:
    file dedupped from ch_for_filtering

    output:
    file "${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta" into (ch_fasta_for_igblast,ch_fasta_for_igblast_filter)

    script:
    """
    SplitSeq.py group -s $dedupped -f CONSCOUNT --num 2 --outname ${dedupped.baseName}_UMI_R1_R2
    sed -n '1~4s/^@/>/p;2~4p' ${dedupped.baseName}_UMI_R1_R2_atleast-2.fastq > ${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta
    """
}

//Run IGBlast
process igblast{
    tag "${fasta.baseName}"

    input:
    file fasta name 'input_igblast.fasta' from ch_fasta_for_igblast
    file igblast from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix) 

    output:
    file "*igblast.fmt7" into ch_igblast_filter

    script:
    """
    AssignGenes.py igblast -s $fasta -b $igblast --organism human --loci ig --format blast
    """
}

//Process output of IGBLAST, makedb + remove non-functional sequences, filter heavy chain and export records to FastA
process igblast_filter {
    tag "${blast.baseName}"
    

    input: 
    file blast name 'blast.fmt7' from ch_igblast_filter
    file fasta name 'fasta.fasta' from ch_fasta_for_igblast_filter
    file imgtbase from ch_imgt_db_for_igblast_filter.mix(ch_imgt_db_for_igblast_filter_mix)

    output:
    file "${blast.baseName}_parse-select.tab" into ch_for_shazam
    file "${blast.baseName}_parse-select_sequences.fasta"

    script:
    """
    MakeDb.py igblast -i ${blast} -s ${fasta} -r ${imgtbase}/human/vdj/imgt_human_IGHV.fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta --regions --scores
    ParseDb.py split -d ${blast.baseName}_db-pass.tab -f FUNCTIONAL
    ParseDb.py select -d ${blast.baseName}_db-pass_FUNCTIONAL-T.tab -f V_CALL -u IGHV --regex --outname ${blast.baseName}
    ConvertDb.py fasta -d ${blast.baseName}_parse-select.tab --if SEQUENCE_ID --sf SEQUENCE_IMGT --mf V_CALL DUPCOUNT
    """
}



/*
    
*/

//Shazam! 
process shazam{
    tag "${tab.baseName}"    
    publishDir "${params.outdir}/shazam/$id", mode: 'copy'

    input:
    file tab from ch_for_shazam
    file imgtbase from ch_imgt_db_for_shazam.mix(ch_imgt_db_for_shazam_mix)
    val id from ch_sample_for_shazam

    output:
    file "threshold.txt" into ch_threshold_for_clone_definition
    file "igh_genotyped.tab" into ch_genotyped_tab_for_clone_definition
    file "v_genotype.fasta" into ch_genotype_fasta_for_germline
    file "Hamming_distance_threshold.pdf" 
    file "genotype.pdf"

    script:
    """
    TIgGER-shazam.R $tab ${imgtbase}/human/vdj/imgt_human_IGHV.fasta
    """
}

//Assign clones
process assign_clones{
    tag "${geno.baseName}" 

    input:
    file geno from ch_genotyped_tab_for_clone_definition
    val threshold from ch_threshold_for_clone_definition

    output:
    file "${geno.baseName}_clone-pass.tab" into ch_for_germlines

    script:
    thr = file(threshold).text
    thr = thr.trim()
    """
    DefineClones.py -d $geno --act set --model ham --norm len --dist $thr --outname ${geno.baseName}
    """
}

//Reconstruct germline sequences
process germline_sequences{
    tag "${clones.baseName}"

    input: 
    file clones from ch_for_germlines
    file imgtbase from ch_imgt_db_for_germline_sequences.mix(ch_imgt_db_for_germline_sequences_mix)
    file geno_fasta from ch_genotype_fasta_for_germline

    output:
    file "${clones.baseName}_germ-pass.tab" into ch_for_alakazam

    script:
    """
    CreateGermlines.py -d ${clones} -g dmask --cloned -r $geno_fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta
    """
}

//Alakazam!
process alakazam{
    tag "${tab.baseName}"
    publishDir "${params.outdir}/alakazam/$id", mode: 'copy'


    input:
    file tab from ch_for_alakazam
    val id from ch_sample_for_alakazam

    output:
    file "*.pdf"
    file "$tab"

    script:
    """
    alakazam.R $tab
    """
}

//Useful functions

 // Return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}