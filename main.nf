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
    log.info nfcoreHeader()
    log.info"""

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

    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
      --imgtdb_base                 Path to predownloaded imgtDB 
      --igblast_base                Path to predownloaded igblastDB

    Define clones:
      --set_cluster_threshold       Set this parameter to allow manual hamming distance threshold for cell cluster definition.
      --cluster_threshold           Once set_cluster_threshold is true, set cluster_threshold value (float).
      --define_clones_only          If set, expects tables as produced by change-O as input, tab delimited. Only performs clonal and germline assignment.
      --changeo_tsv                 Once define_clones_only is specified, set the path to the change-O tsv file.

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
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
   // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// If paths to DBS are provided 
if( params.igblast_base ){
    Channel.fromPath("${params.igblast_base}")
    .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
    .set { ch_igblast_db_for_process_igblast_mix }
} else {
    ch_igblast_db_for_process_igblast_mix = Channel.empty()
}
if( params.imgtdb_base ){
    Channel.fromPath("${params.imgtdb_base}")
    .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
    .into { ch_imgt_db_for_igblast_filter_mix;ch_imgt_db_for_shazam_mix;ch_imgt_db_for_germline_sequences_mix }
} else {
    ch_imgt_db_for_igblast_filter_mix = Channel.empty()
    ch_imgt_db_for_germline_sequences_mix = Channel.empty()
    ch_imgt_db_for_shazam_mix = Channel.empty()
}

//Other parameters
filterseq_q = 20

//Cluster threshold settings
if (params.set_cluster_threshold){
    params.cluster_threshold = 0.0
}

//Define clones only

if (params.define_clones_only){
    params.changeo_tsv = params.changeo_tsv ?: { log.error "No changeo data provided. Make sure you have used the '--changeo_tsv' option."; exit 1 }()
    Channel
        .fromPath( params.changeo_tsv )
        .ifEmpty { exit 1, "Cannot find any changeo_tsv matching: ${params.changeo_tsv}\nNB: Path needs to be enclosed in quotes!" }
        .map { tsv_file -> [file(tsv_file), "${tsv_file.baseName}"] }
        .into { ch_input_tsvs }
} else {
    ch_input_tsvs = Channel.empty()
}

//Set up channels for input primers

if  (!params.define_clones_only){
    Channel.fromPath("${params.cprimers}")
           .ifEmpty{exit 1, "Please specify CPRimers FastA File!"}
           .set {ch_cprimers_fasta}
    Channel.fromPath("${params.vprimers}")
           .ifEmpty{exit 1, "Please specify VPrimers FastA File!"}
           .set { ch_vprimers_fasta }
} else {
    ch_cprimers_fasta = Channel.empty()
    ch_vprimers_fasta = Channel.empty()
}

/*
 * Create a channel for metadata and raw files
 * Columns = id, source, treatment, extraction_time, population, R1, R2, I1
 */
 if  (!params.define_clones_only){
     file_meta = file(params.metadata)
     ch_read_files_for_merge_r1_umi = Channel.from(file_meta)
                    .splitCsv(header: true, sep:'\t')
                    .map { col -> tuple("${col.ID}", "${col.Source}", "${col.Treatment}","${col.Extraction_time}","${col.Population}",returnFile("${col.R1}"),returnFile("${col.R2}"),returnFile("${col.I1}"))}
                    .dump()
 } else {
     ch_read_files_for_merge_r1_umi = Channel.empty()
 }

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision

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
summary['Define clones only']  = params.define_clones_only
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Script dir']     = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-bcellmagic-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/bcellmagic Workflow Summary'
    section_href: 'https://github.com/nf-core/bcellmagic'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

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
    echo "Fetching databases..."

    wget https://raw.githubusercontent.com/nf-core/test-datasets/bcellmagic/database-cache/databases.zip

    unzip databases.zip

    echo "FetchDBs process finished."
    """
}

//Merge I1 UMIs into R1 file
process merge_r1_umi {
    tag "${id}"

    when:
    !params.define_clones_only

    input:
    set val(id), val(source), val(treatment), val(extraction_time), val(population), file(R1), file(R2), file(I1) from ch_read_files_for_merge_r1_umi

    output:
    set file("*UMI_R1.fastq"), file("${R2.baseName}"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_fastqs_for_processing_umi

    script:
    """
    merge_R1_umi.py -R1 "${R1}" -I1 "${I1}" -o "${R1.baseName}_UMI_R1.fastq.gz"
    gunzip "${R1.baseName}_UMI_R1.fastq.gz"
    gunzip -f "${R2}"
    """
}


//Filter by Sequence Quality
process filter_by_sequence_quality {
    tag "${id}"
    publishDir "${params.outdir}/filter_by_sequence_quality/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_fastqs_for_processing_umi

    output:
    set file("${umi.baseName}_quality-pass.fastq"), file("${r2.baseName}_quality-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_filtered_by_seq_quality_for_primer_Masking_UMI
    file "${umi.baseName}_UMI_R1.log"
    file "${r2.baseName}_R2.log"
    file "${umi.baseName}_UMI_R1_table.tab"
    file "${r2.baseName}_R2_table.tab"
    file "command_log.txt"

    script:
    """
    FilterSeq.py quality -s $umi -q $filterseq_q --outname "${umi.baseName}" --log "${umi.baseName}_UMI_R1.log"
    FilterSeq.py quality -s $r2 -q $filterseq_q --outname "${r2.baseName}" --log "${r2.baseName}_R2.log"
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${umi.baseName}_UMI_R1.log" "${r2.baseName}_R2.log" -f ID QUALITY
    """
}

//Mask them primers
process mask_primers {
    tag "${id}"
    publishDir "${params.outdir}/mask_primers/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi_file), file(r2_file), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_filtered_by_seq_quality_for_primer_Masking_UMI
    file(cprimers) from ch_cprimers_fasta.collect() 
    file(vprimers) from ch_vprimers_fasta.collect()

    output:
    set file("${umi_file.baseName}_UMI_R1_primers-pass.fastq"), file("${r2_file.baseName}_R2_primers-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_pair_seq_umi_file
    file "${umi_file.baseName}_UMI_R1.log"
    file "${r2_file.baseName}_R2.log"
    file "command_log.txt"

    script:
    """
    MaskPrimers.py score --nproc ${task.cpus} -s $umi_file -p ${cprimers} --start 8 --mode cut --barcode --outname ${umi_file.baseName}_UMI_R1 --log ${umi_file.baseName}_UMI_R1.log
    MaskPrimers.py score --nproc ${task.cpus} -s $r2_file -p ${vprimers} --start 0 --mode mask --outname ${r2_file.baseName}_R2 --log ${r2_file.baseName}_R2.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${umi_file.baseName}_UMI_R1.log" "${r2_file.baseName}_R2.log" -f ID PRIMER ERROR
    """
}

//Pair the UMI_R1 and R2
process pair_seq{
    tag "${id}"
    publishDir "${params.outdir}/pair_sequences/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_pair_seq_umi_file

    output:
    set file("${umi.baseName}_pair-pass.fastq"), file("${r2.baseName}_pair-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_umi_for_umi_cluster_sets
    file "command_log.txt"

    script:
    """
    PairSeq.py -1 $umi -2 $r2 --1f BARCODE --coord illumina
    cp ".command.out" "command_log.txt"
    """
}

//Cluster sets to deal with too low UMI diversity
process cluster_sets {
    tag "${id}"
    publishDir "${params.outdir}/cluster_sets/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_umi_for_umi_cluster_sets

    when:
    !params.define_clones_only

    output:
    set file ("${umi.baseName}_UMI_R1_cluster-pass.fastq"), file("${r2.baseName}_R2_cluster-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_umi_for_reheader
    file "command_log.txt"

    script:
    """
    ClusterSets.py set --nproc ${task.cpus} -s $umi --outname ${umi.baseName}_UMI_R1 
    ClusterSets.py set --nproc ${task.cpus} -s $r2 --outname ${r2.baseName}_R2
    cp ".command.out" "command_log.txt"
    """
}

//ParseHeaders to annotate barcode into cluster field
process reheader {
    tag "${id}"

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_umi_for_reheader

    output:
    set file("${umi.baseName}_reheader.fastq"), file("${r2.baseName}_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_umi_for_consensus

    script:
    """
    ParseHeaders.py copy -s $umi -f BARCODE -k CLUSTER --act cat
    ParseHeaders.py copy -s $r2 -f BARCODE -k CLUSTER --act cat 
    """
}


//Build UMI consensus
process build_consensus{
    tag "${id}"
    publishDir "${params.outdir}/build_consensus/$id", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_umi_for_consensus

    output:
    set file("${umi.baseName}_UMI_R1_consensus-pass.fastq"), file("${r2.baseName}_R2_consensus-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_consensus_passed_umi
    file "${umi.baseName}_UMI_R1.log"
    file "${r2.baseName}_R2.log"
    file "${umi.baseName}_UMI_R1_table.tab"
    file "${r2.baseName}_R2_table.tab"
    file "command_log.txt"

    script:
    """
    BuildConsensus.py -s $umi --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${umi.baseName}_UMI_R1 --log ${umi.baseName}_UMI_R1.log
    BuildConsensus.py -s $r2 --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${r2.baseName}_R2 --log ${r2.baseName}_R2.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${umi.baseName}_UMI_R1.log" "${r2.baseName}_R2.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
    """
}

//Re-pair UMI_R1+R2
process repair{
    tag "${id}"
    publishDir "${params.outdir}/repair_mates/$id", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_consensus_passed_umi

    output:
    set file("*UMI_R1_consensus-pass_pair-pass.fastq"), file("*R2_consensus-pass_pair-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_repaired_UMI_for_assembly
    file "command_log.txt"

    script:
    """
    PairSeq.py -1 $umi -2 $r2 --coord presto
    cp ".command.out" "command_log.txt"
    """
}
   

//Assemble the mate pairs
process assemble{
    tag "${id}"
    publishDir "${params.outdir}/assemble_pairs/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(umi), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_repaired_UMI_for_assembly

    output:
    set file("${umi.baseName}_UMI_R1_R2_assemble-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_combine_UMI
    file "${umi.baseName}_UMI_R1_R2.log"
    file "${umi.baseName}_UMI_R1_R2_table.tab"
    file "command_log.txt"

    script:
    """
    AssemblePairs.py align -1 $umi -2 $r2 --coord presto --rc tail --1f CONSCOUNT PRCONS --2f CONSCOUNT PRCONS --outname ${umi.baseName}_UMI_R1_R2 --log ${umi.baseName}_UMI_R1_R2.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${umi.baseName}_UMI_R1_R2.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT LENGTH OVERLAP ERROR PVALUE
    """
}

//combine UMI read group size annotations
process combine_umi_read_groups{
    tag "${id}" 

    when:
    !params.define_clones_only

    input:
    set file(assembled), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_combine_UMI

    output:
    set file("*_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_prcons_parseheaders

    script:
    """
    ParseHeaders.py collapse -s $assembled -f CONSCOUNT --act min
    """
}


//Copy field PRCONS to have annotation for C_primer and V_primer independently
process copy_prcons{
    tag "${id}" 

    when:
    !params.define_clones_only

    input:
    set file(combined), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_prcons_parseheaders

    output:
    set file("*reheader_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_metadata_anno

    script:
    """
    ParseHeaders.py copy -s $combined -f PRCONS PRCONS --act first last -k C_PRIMER V_PRIMER
    """
}

//Add Metadata annotation to headers
process metadata_anno{
    tag "${id}"

    when:
    !params.define_clones_only

    input:
    set file(prcons), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_metadata_anno

    output:
    set file("*_reheader_reheader_reheader.fastq"), val("$id"), val("$source") into ch_for_dedup

    script:
    """
    ParseHeaders.py add -s $prcons -f SAMPLE_CODE SOURCE TREATMENT EXTRACT_TIME POPULATION -u $id $source $treatment $extraction_time $population
    """
}

//Removal of duplicate sequences
process dedup {
    tag "${id}"
    publishDir "${params.outdir}/deduplicates/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "fastq/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(dedup), val(id), val(source) from ch_for_dedup

    output:
    set file("${dedup.baseName}_UMI_R1_R2_collapse-unique.fastq"), val("$id"), val("$source") into ch_for_filtering
    file "${dedup.baseName}_UMI_R1_R2.log"
    file "${dedup.baseName}_UMI_R1_R2_table.tab"
    file "command_log.txt"

    script:
    """
    CollapseSeq.py -s $dedup -n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname ${dedup.baseName}_UMI_R1_R2 --log ${dedup.baseName}_UMI_R1_R2.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${dedup.baseName}_UMI_R1_R2.log" -f HEADER DUPCOUNT
    """
}

//Filtering to sequences with at least two representative reads and convert to FastA
process filter_seqs{
    tag "${id}"
    publishDir "${params.outdir}/filter_representative_2/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fasta") > 0) "fasta/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input:
    set file(dedupped), val(id), val(source) from ch_for_filtering

    output:
    set file("${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta"), val("$id"), val("$source") into ch_fasta_for_igblast
    file "${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta"
    file "command_log.txt"

    script:
    """
    SplitSeq.py group -s $dedupped -f CONSCOUNT --num 2 --outname ${dedupped.baseName}_UMI_R1_R2
    cp ".command.out" "command_log.txt"
    sed -n '1~4s/^@/>/p;2~4p' ${dedupped.baseName}_UMI_R1_R2_atleast-2.fastq > ${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta
    """
}

//Run IGBlast
process igblast{
    tag "${id}"

    when:
    !params.define_clones_only

    input:
    set file('input_igblast.fasta'), val(id), val(source) from ch_fasta_for_igblast
    file igblast from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect() 

    output:
    set file("*igblast.fmt7"), file('input_igblast.fasta'), val("$id"), val("$source") into ch_igblast_filter

    script:
    """
    AssignGenes.py igblast -s input_igblast.fasta -b $igblast --organism human --loci ig --format blast
    """
}


//Process output of IGBLAST, makedb + remove non-functional sequences, filter heavy chain and export records to FastA
process igblast_filter {
    tag "${id}"
    publishDir "${params.outdir}/igblast/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fasta") > 0) "fasta/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename.indexOf(".tab") > 0) "table/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    when:
    !params.define_clones_only

    input: 
    set file('blast.fmt7'), file('fasta.fasta'), val(id), val(source) from ch_igblast_filter
    file imgtbase from ch_imgt_db_for_igblast_filter.mix(ch_imgt_db_for_igblast_filter_mix).collect()

    output:
    set file("${base}_parse-select.tab"), val("$id"), val("$source") into ch_for_merge
    file "${base}_parse-select_sequences.fasta"
    file "${base}_parse-select.tab"
    file "command_log.txt"

    script:
    base = "blast"
    """
    MakeDb.py igblast -i blast.fmt7 -s fasta.fasta -r ${imgtbase}/human/vdj/imgt_human_IGHV.fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta --regions --scores
    ParseDb.py split -d ${base}_db-pass.tab -f FUNCTIONAL
    ParseDb.py select -d ${base}_db-pass_FUNCTIONAL-T.tab -f V_CALL -u IGHV --regex --outname ${base}
    ConvertDb.py fasta -d ${base}_parse-select.tab --if SEQUENCE_ID --sf SEQUENCE_IMGT --mf V_CALL DUPCOUNT
    cp ".command.out" "command_log.txt"
    """
}

//Merge tables belonging to the same patient
process merge_tables{
    tag "merge tables"
    publishDir "${params.outdir}/shazam/$source", mode: 'copy'

    input:
    set file(tab), val(id), val(source) from ch_for_merge.groupTuple(by: 3)

    output:
    set file("tab"), val(source) into ch_for_shazam

    script:
    """
    """

}


//Shazam! 
process shazam{
    tag "${id}"    
    publishDir "${params.outdir}/shazam/$id", mode: 'copy'

    input:
    set file(tab), val(id) from ch_for_shazam.mix(ch_input_tsvs)
    file imgtbase from ch_imgt_db_for_shazam.mix(ch_imgt_db_for_shazam_mix).collect()

    output:
    set file("threshold.txt"), file("igh_genotyped.tab"), file("v_genotype.fasta"), val("$id") into ch_threshold_for_clone_definition
    file "Hamming_distance_threshold.pdf" 
    file "genotype.pdf"

    script:
    """
    TIgGER-shazam.R $tab ${imgtbase}/human/vdj/imgt_human_IGHV.fasta
    """
}

//Assign clones
process assign_clones{
    tag "${id}" 
    publishDir "${params.outdir}/define_clones/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fasta") > 0) "fasta/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename.indexOf(".tab") > 0) "table/$filename"
            else if (filename == "command_log.txt") "$filename"
            else if (filename == "threshold.txt" && !params.set_cluster_threshold) "$filename"
            else null
        }

    input:
    set val(threshold), file(geno), file(geno_fasta), val(id) from ch_threshold_for_clone_definition

    output:
    set file("${geno.baseName}_clone-pass.tab"), file("$geno_fasta"), val("$id") into ch_for_germlines
    file "${geno.baseName}_clone-pass.tab"
    file "${geno.baseName}_table.tab"
    file "command_log.txt"

    script:
    if (params.set_cluster_threshold) {
        thr = params.cluster_threshold
    } else {
        thr = file(threshold).text
        thr = thr.trim()
    }
    """
    DefineClones.py -d $geno --act set --model ham --norm len --dist $thr --outname ${geno.baseName} --log ${geno.baseName}.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${geno.baseName}.log" -f ID VCALL JCALL JUNCLEN CLONED FILTERED CLONES
    """
}


//Reconstruct germline sequences
process germline_sequences{
    tag "${id}"
    publishDir "${params.outdir}/create_germlines/$id", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fasta") > 0) "fasta/$filename"
            else if (filename.indexOf(".log") > 0) null
            else if (filename.indexOf("table.tab") > 0) "info/$filename"
            else if (filename.indexOf(".tab") > 0) "table/$filename"
            else if (filename == "command_log.txt") "$filename"
            else null
        }

    input: 
    set file(clones), file(geno_fasta), val(id) from ch_for_germlines
    file imgtbase from ch_imgt_db_for_germline_sequences.mix(ch_imgt_db_for_germline_sequences_mix).collect()

    output:
    set file("${clones.baseName}_germ-pass.tab"), val("$id") into ch_for_alakazam
    file "${clones.baseName}_germ-pass.tab"
    file "command_log.txt"

    script:
    """
    CreateGermlines.py -d ${clones} -g dmask --cloned -r $geno_fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta --log ${clones.baseName}.log
    cp ".command.out" "command_log.txt"
    ParseLog.py -l "${clones.baseName}.log" -f ID V_CALL D_CALL J_CALL
    """
}

//Alakazam!
process alakazam{
    tag "${id}" 
    publishDir "${params.outdir}/alakazam/$id", mode: 'copy'


    input:
    set file(tab), val(id) from ch_for_alakazam

    output:
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

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/bcellmagic] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/bcellmagic] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/bcellmagic] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/bcellmagic] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/bcellmagic] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/bcellmagic] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/bcellmagic]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/bcellmagic]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/bcellmagic v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
