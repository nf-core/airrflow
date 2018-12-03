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
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

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
output_docs = Channel.fromFile("$baseDir/docs/output.md")


//Set up channels for input primers
Channel.fromPath("${params.cprimers}")
       .ifEmpty{exit 1, "Please specify CPRimers FastA File!"}
       .into {ch_cprimers_fasta}
Channel.fromPath("${params.vprimers}")
       .ifEmpty{exit 1, "Please specify VPrimers FastA File!"}
       .into { ch_vprimers_fasta }

//Check for supplied databases using groovy syntax
params.igblast_base ? ch_igblast = Channel.fromPath(igblast_base, checkIfExists: true) : igblast_base = false

params.imgt_base ? ch_imgt_base = Channel.fromPath(imgt_base, checkIfExists: true) : imgt_base = false

saveDBs = false

//Other parameters
filterseq_q=20


//Download data process
process fetchDBs{
    tag "fetchBlastDBs"

    publishDir path: { params.saveDBs ? "${params.outdir}/dbs" : params.outdir },
    saveAs: { params.saveDBs ? it : null }, mode: 'copy'

    when: (!"${params.igblast_base}" || !"${params.imgt_base}")

    output:
    file "$igblast_base" into ch_igblast_for_mixing
    file "$imgt_base" into ch_imgt_for_mixing
    
    script:
    """
    fetch_igblastdb.sh -o $igblast_basew
    fetch_imgtdb.sh -o $imgt_base
    imgt2igblast.sh -i $imgt_base -o $igblast_base
    """
}

//Mix channels for DBs
//Igblast DB
ch_igblast.mix(ch_igblast_for_mixing).into { ch_igblast_db_for_process_igblast }
//IMGT DB
ch_imgt_base.mix(ch_imgt_for_mixing).into {ch_imgt_db_for_igblast_filter;ch_imgt_db_for_shazam;ch_imgt_db_for_germline_sequences}

/*
 * Create a channel for input read files
 */
 if(params.readPaths){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { ch_read_files_for_merge_r1_umi }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, 3 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\n" }
         .into { ch_read_files_for_merge_r1_umi }
 }

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
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
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


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

//Merge I1 UMIs into R1 file
process merge_r1_umi {
    tag $read 

    input:
    set val(name), file(reads) from ch_read_files_for_merge_r1_umi 

    output:
    file "{*_UMI_*,*_R2_*}.fastq.gz" into ch_umi_merged_for_decompression

    script:
    """
    merge_R1_umi.py -R1 $reads[0] -I1 $reads[3] -o "${reads[0].baseName}_UMI_R1.fastq.gz"
    """
}

//Decompress all the stuff
process decompress {
    tag "$reads[0]"
    
    input:
    file(reads) from ch_umi_merged_for_decompression

    output:
    file "*.fastq" into ch_fastqs_for_processing

    script:
    """
    gunzip *.fastq.gz 
    """
}

//Filter by Sequence Quality
process filter_by_sequence_quality {
    tag "$reads[0]"

    input:
    file(reads) fom ch_fastqs_for_processing

    output:
    file "${reads[0]}*quality-pass.fastq" into ch_filtered_by_seq_quality_for_primer_Masking_UMI
    file "${reads[1]}*quality-pass.fastq" into ch_filtered_by_seq_quality_for_primerMasking_R2

    script:
    """
    FilterSeq.py quality -s "$reads[0]" -q $filterseq_q --outname "${reads[0].baseName}"
    FilterSeq.py quality -s "$reads[1]" -q $filterseq_q --outname "${reads[1].baseName}"
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
    BuildConsensus.py -s $umi --bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${umi.baseName}_UMI_R1
    BuildConsensus.py -s $r2 --bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${r2.baseName}_R2
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
//TODO check TREATMENT,EXTRACT_TIME and POPULATION as well! 
process metadata_anno{
    tag "${prcons.baseName}"

    input:
    file prcons from ch_for_metadata_anno

    output:
    file "*_reheader_reheader_reheader.fastq" into ch_for_dedup 

    script:
    """
    ParseHeaders.py add -s $prcons -f TREATMENT EXTRACT_TIME POPULATION -u ${TREATMENT} ${EXTRACT_TIME} ${POPULATION}
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
    file "${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta" into ch_fasta_for_igblast

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
    file fasta from ch_fasta_for_igblast
    file igblast from ch_igblast_db_for_process_igblast 

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
    file blast from ch_igblast_filter
    file imgtbase from ch_imgt_db_for_igblast_filter

    output:
    file "${blast.baseName}_UMI_R1_R2_atleast-2_igblast_db-pass_FUNCTIONAL-T_parse-select.tab" into ch_for_shazam

    script:
    """
    MakeDb.py igblast -i $blast -s ${blast.baseName}_UMI_R1_R2_atleast-2.fasta -r ${imgtbase}/human/vdj/imgt_human_IGHV.fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta --regions --scores
    ParseDb.py split -d ${blast.baseName}_UMI_R1_R2_atleast-2_igblast_db-pass.tab -f FUNCTIONAL
    ParseDb.py select -d ${blast.baseName}_UMI_R1_R2_atleast-2_igblast_db-pass_FUNCTIONAL-T.tab -f V_CALL -u IGHV --regex --outname ${blast.baseName}_UMI_R1_R2_atleast-2_igblast_db-pass_FUNCTIONAL-T
    ConvertDb.py fasta -d ${blast.baseName}_UMI_R1_R2_atleast-2_igblast_db-pass.tab --if SEQUENCE_ID --sf SEQUENCE_IMGT --mf V_CALL DUPCOUNT
    """
}

//Shazam! 
process shazam{
    tag "${tab.baseName}"

    input:
    file tab from ch_for_shazam
    file imgtbase from ch_imgt_db_for_shazam

    output:
    file "threshold.txt" into ch_threshold_for_clone_definition
    file "igh_genotyped.tab" into ch_genotyped_tab_for_clone_definition
    file "v_genotype.fasta" into ch_genotype_fasta_for_germline

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
    file threshold from ch_threshold_for_clone_definition

    output:
    file "${geno.baseName}_clone-pass.tab" into ch_for_germlines

    script:
    """
    thr=`cat $threshold`
    DefineClones.py -d $geno --act set --model ham --norm len --dist ${thr} --outname
    ${geno.baseName}
    """
}

//Reconstruct germline sequences
process germline_sequences{
    tag "${clones.baseName}"

    input: 
    file clones from ch_for_germlines
    file imgtbase from ch_imgt_db_for_germline_sequences
    file geno_fasta from ch_genotype_fasta_for_germline

    output:
    file "${clones.baseName}_clone-pass.tab" into ch_for_alakazam

    script:
    """
    CreateGermlines.py -d ${clones.baseName} -g dmask --cloned -r $geno_fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta
    """
}

//Alakazam!
process alakazam{
    tag "${tab.baseName}"

    input:
    file tab from ch_for_alakazam

    output:
    file "*"

    script:
    """
    alakazam.R $tab
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from output_docs

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
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
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
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/bcellmagic] Pipeline Complete"

}
