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

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    // TODO nf-core: Update typical command used to run pipeline
    def command = "nextflow run nf-core/bcellmagic -profile <docker/singularity/podman/conda/institute> --input metadata.tsv --cprimers CPrimers.fasta --vprimers VPrimers.fasta"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, "Please provide input file containing the sample metadata with the '--input' option." }


// // If paths to DBS are provided 
// if( params.igblast_base ){
//     Channel.fromPath("${params.igblast_base}")
//     .ifEmpty { exit 1, "IGBLAST DB not found: ${params.igblast_base}" }
//     .set { ch_igblast_db_for_process_igblast_mix }
// } else {
//     ch_igblast_db_for_process_igblast_mix = Channel.empty()
// }
// if( params.imgtdb_base ){
//     Channel.fromPath("${params.imgtdb_base}")
//     .ifEmpty { exit 1, "IMGTDB not found: ${params.imgtdb_base}" }
//     .into { ch_imgt_db_for_igblast_filter_mix;ch_imgt_db_for_shazam_mix;ch_imgt_db_for_germline_sequences_mix }
// } else {
//     ch_imgt_db_for_igblast_filter_mix = Channel.empty()
//     ch_imgt_db_for_germline_sequences_mix = Channel.empty()
//     ch_imgt_db_for_shazam_mix = Channel.empty()
// }


// Validate primer protocol
if (params.cprimers)  { ch_cprimers_fasta = Channel.fromPath(params.cprimers, checkIfExists: true) } else { exit 1, "Please provide a C-region primers fasta file with the '--cprimers' option." }
if (!params.race_5prime){
    if (params.vprimers)  { 
        ch_vprimers_fasta = Channel.fromPath(params.vprimers, checkIfExists: true) 
    } else { 
        exit 1, "Please provide a V-region primers fasta file with the '--vprimers' option, or specify a 5'RACE protocol with the '--race_5prime' option." 
    }
} else {
    if (params.vprimers) { 
        exit 1, "The 5' RACE protocol does not accept V-region primers, please remove the option '--vprimers' or the option '--race_5prime'."
    } else if (params.race_linker) {
        ch_vprimers_fasta = Channel.fromPath(params.race_linker, checkIfExists: true)
    } else {
        exit 1, "The 5' RACE protocol requires a linker or Template Switch oligo sequence, please provide it with the option '--race_linker'."
    }
}

// //Validate UMI position
// if (params.index_file & params.umi_position == 'R2') {exit 1, "Please do not set `--umi_position` option if index file with UMIs is provided."}
// if (params.umi_length == 0) {exit 1, "Please provide the UMI barcode length in the option `--umi_length`."}
// if (!params.index_file & params.umi_start != 0) {exit 1, "Setting a UMI start position is only allowed when providing the UMIs in a separate index read file. If so, please provide the `--index_file` flag as well."}
// //Read processed tabs if downstream_only 
// if (params.downstream_only){
//     Channel
//         .fromFilePairs(params.changeo_tables, size: 1) {file -> file.baseName}
//         .ifEmpty {exit 1, "Cannot find any changeo tables matching: ${params.changeo_tables}.\nTry enclosing paths in quotes!\nTry adding a * wildcard!"}
//         .set {ch_tabs_for_clonal_analysis}
//         .println()
// } else {
//     ch_tabs_for_clonal_analysis = Channel.empty()
// }

/*
 * Create a channel for metadata and raw files
 * Columns = id, source, treatment, extraction_time, population, R1, R2, I1
 */

// if (!params.downstream_only){
//     Channel.fromPath("${params.input}")
//             .ifEmpty{exit 1, "Please provide metadata file!"}
//             .set { ch_metadata_file_for_process_logs }

//     if (params.index_file) {
//         ch_read_files_for_merge_r1_umi_index = Channel.from( ch_metadata )
//                 .splitCsv(header: true, sep:'\t')
//                 .map { col -> tuple("${col.ID}", "${col.Source}", "${col.Treatment}","${col.Extraction_time}","${col.Population}", file("${col.R1}", checkifExists: true),file("${col.R2}", checkifExists: true), file("${col.I1}", checkifExists: true))}
//                 .dump()
//         ch_read_files_for_merge_r1_umi = Channel.empty()
//     } else {
//         ch_read_files_for_merge_r1_umi = Channel.from( ch_metadata )
//                 .splitCsv(header: true, sep:'\t')
//                 .map { col -> tuple("${col.ID}", "${col.Source}", "${col.Treatment}","${col.Extraction_time}","${col.Population}", file("${col.R1}", checkifExists: true), file("${col.R2}", checkifExists: true))}
//                 .dump()
//         ch_read_files_for_merge_r1_umi_index = Channel.empty()
//     }
// } else {
//     ch_read_files_for_merge_r1_umi = Channel.empty()
//     ch_read_files_for_merge_r1_umi_index = Channel.empty()
//     ch_metadata_file_for_process_logs = Channel.empty()
// }


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions'  addParams( options: [publish_files : ['csv':'']] )
include { MERGE_UMI } from './modules/local/process/merge_UMI'                          addParams( options: [:] )
include { GUNZIP } from './modules/local/process/gunzip'                                addParams( options: [:] )
include { PRESTO_FILTERSEQ } from './modules/local/process/presto_filterseq'            addParams( options: modules['presto_filterseq'] )
include { PRESTO_MASKPRIMERS } from './modules/local/process/presto_maskprimers'        addParams( options: modules['presto_maskprimers'] )
include { PRESTO_PAIRSEQ } from './modules/local/process/presto_pairseq'                addParams( options: modules['presto_pairseq'] )
include { PRESTO_CLUSTERSETS } from './modules/local/process/presto_clustersets'        addParams( options: modules['presto_clustersets'] )

// Local: Sub-workflows
include { INPUT_CHECK           } from './modules/local/subworkflow/input_check'       addParams( options: [:] )

// nf-core/modules: Modules
include { FASTQC                } from './modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc'] )
include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input )
    .groupTuple(by: [0])
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_fastqc }

    ch_merge_umi_gunzip = ch_fastqc.map{ it -> it.flatten() }

    //FastQC
    FASTQC (ch_fastqc)
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    
    //Merge UMI from index file to R1 if provided
    if (params.index_file) {
        MERGE_UMI ( ch_merge_umi_gunzip )
        .set{ ch_gunzip }
    } else {
        ch_gunzip = ch_merge_umi_gunzip
    }

    //GUNZIP: gunzip fastq.gz to fastq
    GUNZIP ( ch_gunzip )
    ch_software_versions = ch_software_versions.mix(GUNZIP.out.version.first().ifEmpty(null))

    //PRESTO FILTERSEQ: Filter sequences by quality score
    PRESTO_FILTERSEQ ( GUNZIP.out.reads )
    ch_software_versions = ch_software_versions.mix(PRESTO_FILTERSEQ.out.version.first().ifEmpty(null))

    //PRESTO MASKPRIMERS: Mask primers
    PRESTO_MASKPRIMERS ( 
        PRESTO_FILTERSEQ.out.reads,
        ch_cprimers_fasta.collect(),
        ch_vprimers_fasta.collect()
    )

    //PRESTO PAIR: Pre-consensus pair
    PRESTO_PAIRSEQ (
        PRESTO_MASKPRIMERS.out.reads
    )

    //PRESTO CLUSTERSETS: cluster sequences by similarity
    PRESTO_CLUSTERSETS (
        PRESTO_PAIRSEQ.out.reads
    )
    ch_software_versions = ch_software_versions.mix(PRESTO_CLUSTERSETS.out.version.first().ifEmpty(null))

    // Software versions
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        
        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////



//////////////////////
// Old part of main.nf
//////////////////////


// /*
//  * Parse software version numbers
//  */
// process get_software_versions {
//     publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
//     saveAs: {filename ->
//         if (filename.indexOf(".csv") > 0) filename
//         else null
//     }

//     output:
//     file 'software_versions_mqc.yaml' into ch_software_versions_yaml
//     file "software_versions.csv"

//     script:
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version &> v_fastqc.txt
//     multiqc --version &> v_multiqc.txt
//     vsearch --version &> v_vsearch.txt
//     muscle -version &> v_muscle.txt
//     python -c "import presto; print(presto.__version__)" > v_presto.txt
//     python -c "import changeo; print(changeo.__version__)" > v_changeo.txt
//     echo \$(R --version 2>&1) > v_R.txt
//     Rscript -e "library(shazam); write(x=as.character(packageVersion('shazam')), file='v_shazam.txt')"
//     Rscript -e "library(alakazam); write(x=as.character(packageVersion('alakazam')), file='v_alakazam.txt')"
//     Rscript -e "library(tigger); write(x=as.character(packageVersion('tigger')), file='v_tigger.txt')"
//     scrape_software_versions.py &> software_versions_mqc.yaml
//     """
// }

// //Download data process
// process fetchDBs{
//     tag "fetchBlastDBs"

//     publishDir path: { params.save_databases ? "${params.outdir}/dbs" : params.outdir },
//     saveAs: { params.save_databases ? it : null }, mode: params.publish_dir_mode

//     when:
//     !params.igblast_base | !params.imgtdb_base

//     output:
//     file "igblast_base" into ch_igblast_db_for_process_igblast
//     file "imgtdb_base" into (ch_imgt_db_for_igblast_filter,ch_imgt_db_for_shazam,ch_imgt_db_for_germline_sequences)
    
//     script:
//     """
//     echo "Fetching databases..."

//     wget https://raw.githubusercontent.com/nf-core/test-datasets/bcellmagic/database-cache/databases.zip

//     unzip databases.zip

//     echo "FetchDBs process finished."
//     """
// }


// //Cluster sets to deal with too low UMI diversity
// process cluster_sets {
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/cluster_sets/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "$filename"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(r1), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_cluster_sets

//     output:
//     set file ("${r1.baseName}_R1_cluster-pass.fastq"), file("${r2.baseName}_R2_cluster-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_reheader
//     file "${id}_command_log.txt" into cluster_sets_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     ClusterSets.py set --nproc ${task.cpus} -s $r1 --outname ${r1.baseName}_R1 > "${id}_command_log.txt"
//     ClusterSets.py set --nproc ${task.cpus} -s $r2 --outname ${r2.baseName}_R2 >> "${id}_command_log.txt"
//     """
// }

// //ParseHeaders to annotate barcode into cluster field
// process reheader {
//     tag "${id}"

//     input:
//     set file(r1), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_reheader

//     output:
//     set file("${r1.baseName}_reheader.fastq"), file("${r2.baseName}_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_consensus

//     when:
//     !params.downstream_only

//     script:
//     """
//     ParseHeaders.py copy -s $r1 -f BARCODE -k CLUSTER --act cat > "${id}_command_log.txt"
//     ParseHeaders.py copy -s $r2 -f BARCODE -k CLUSTER --act cat >> "${id}_command_log.txt"
//     """
// }


// //Build UMI consensus
// process build_consensus{
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/build_consensus/$id", mode: params.publish_dir_mode,
//     saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "$filename"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(r1), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_consensus

//     output:
//     set file("${r1.baseName}_R1_consensus-pass.fastq"), file("${r2.baseName}_R2_consensus-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_post_consensus_pair
//     file "${r1.baseName}_R1.log"
//     file "${r2.baseName}_R2.log"
//     file "${r1.baseName}_R1_table.tab"
//     file "${r2.baseName}_R2_table.tab"
//     file "${id}_command_log.txt" into build_consensus_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     BuildConsensus.py -s $r1 --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror 0.1 --maxgap 0.5 --outname ${r1.baseName}_R1 --log ${r1.baseName}_R1.log > "${id}_command_log.txt"
//     BuildConsensus.py -s $r2 --bf CLUSTER --nproc ${task.cpus} --pf PRIMER --prcons $params.primer_consensus --maxerror 0.1 --maxgap 0.5 --outname ${r2.baseName}_R2 --log ${r2.baseName}_R2.log >> "${id}_command_log.txt"
//     ParseLog.py -l "${r1.baseName}_R1.log" "${r2.baseName}_R2.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
//     """
// }

// //Re-pair R1 and R2
// process post_consensus_pair{
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/pair_mates_post_consensus/$id", mode: params.publish_dir_mode,
//     saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "${id}_repair_mates_logs.tab"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(r1), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_post_consensus_pair

//     output:
//     set file("*R1_consensus-pass_pair-pass.fastq"), file("*R2_consensus-pass_pair-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_assembly
//     file "${id}_command_log.txt" into repair_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     PairSeq.py -1 $r1 -2 $r2 --coord presto > "${id}_command_log.txt"
//     """
// }
   

// //Assemble the mate pairs
// process assemble_pairs{
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/assemble_pairs/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "${id}_assemble_pairs_logs.tab"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(r1), file(r2), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_assembly

//     output:
//     set file("${id}_assemble-pass.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_combine_UMI
//     file "${id}.log"
//     file "${id}_table.tab"
//     file "${id}_command_log.txt" into assemble_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     AssemblePairs.py align -1 $r1 -2 $r2 --coord presto --rc tail --1f CONSCOUNT PRCONS --2f CONSCOUNT PRCONS --outname ${id} --log ${id}.log > "${id}_command_log.txt"
//     ParseLog.py -l "${id}.log" -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT LENGTH OVERLAP ERROR PVALUE
//     """
// }

// //combine UMI read group size annotations
// process combine_umi_read_groups{
//     tag "${id}"

//     input:
//     set file(assembled), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_combine_UMI

//     output:
//     set file("*_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_prcons_parseheaders

//     when:
//     !params.downstream_only
    
//     script:
//     """
//     ParseHeaders.py collapse -s $assembled -f CONSCOUNT --act min
//     """
// }


// //Copy field PRCONS to have annotation for C_primer and V_primer independently
// process copy_prcons{
//     tag "${id}"

//     input:
//     set file(combined), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_prcons_parseheaders

//     output:
//     set file("*reheader_reheader.fastq"), val("$id"), val("$source"), val("$treatment"), val("$extraction_time"), val("$population") into ch_for_metadata_anno

//     when:
//     !params.downstream_only

//     script:
//     """
//     ParseHeaders.py copy -s $combined -f PRCONS PRCONS --act first last -k C_PRIMER V_PRIMER
//     """
// }

// //Add Metadata annotation to headers
// process metadata_anno{
//     tag "${id}"

//     input:
//     set file(prcons), val(id), val(source), val(treatment), val(extraction_time), val(population) from ch_for_metadata_anno

//     output:
//     set file("*_reheader_reheader_reheader.fastq"), val("$id"), val("$source") into ch_for_dedup

//     when:
//     !params.downstream_only

//     script:
//     """
//     ParseHeaders.py add -s $prcons -f SAMPLE_CODE SOURCE TREATMENT EXTRACT_TIME POPULATION -u $id $source $treatment $extraction_time $population
//     """
// }

// //Removal of duplicate sequences
// process deduplicate {
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/deduplicate/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "${id}_deduplicate_logs.tab"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(dedup), val(id), val(source) from ch_for_dedup

//     output:
//     set file("${dedup.baseName}_UMI_R1_R2_collapse-unique.fastq"), val("$id"), val("$source") into ch_for_filtering
//     file "${dedup.baseName}_UMI_R1_R2.log"
//     file "${dedup.baseName}_UMI_R1_R2_table.tab"
//     file "${id}_command_log.txt" into dedup_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     CollapseSeq.py -s $dedup -n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname ${dedup.baseName}_UMI_R1_R2 --log ${dedup.baseName}_UMI_R1_R2.log > "${id}_command_log.txt"
//     ParseLog.py -l "${dedup.baseName}_UMI_R1_R2.log" -f HEADER DUPCOUNT
//     """
// }

// //Filtering to sequences with at least two representative reads and convert to FastA
// process filter_representative_2{
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/filter_representative_2/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "$filename"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input:
//     set file(dedupped), val(id), val(source) from ch_for_filtering

//     output:
//     set file("${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta"), val("$id"), val("$source") into ch_fasta_for_igblast
//     file "${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta"
//     file "${id}_command_log.txt" into filter_seqs_log

//     when:
//     !params.downstream_only

//     script:
//     """
//     SplitSeq.py group -s $dedupped -f CONSCOUNT --num 2 --outname ${dedupped.baseName}_UMI_R1_R2 > "${id}_command_log.txt"
//     sed -n '1~4s/^@/>/p;2~4p' ${dedupped.baseName}_UMI_R1_R2_atleast-2.fastq > ${dedupped.baseName}_UMI_R1_R2_atleast-2.fasta
//     """
// }

// //Run IGBlast
// process igblast{
//     tag "${id}"

//     input:
//     set file('input_igblast.fasta'), val(id), val(source) from ch_fasta_for_igblast
//     file igblast from ch_igblast_db_for_process_igblast.mix(ch_igblast_db_for_process_igblast_mix).collect() 

//     output:
//     set file("*igblast.fmt7"), file('input_igblast.fasta'), val("$id"), val("$source") into ch_igblast_filter

//     when:
//     !params.downstream_only

//     script:
//     """
//     AssignGenes.py igblast -s input_igblast.fasta -b $igblast --organism $params.species --loci $params.loci --format blast
//     """
// }


// //Process output of IGBLAST, makedb + remove non-functional sequences, filter heavy chain and export records to FastA
// process igblast_filter {
//     tag "${id}"
//     publishDir "${params.outdir}/preprocessing/igblast/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "$filename"
//             else if (filename.indexOf(".fasta") > 0) "fasta/$filename"
//             else if (filename.indexOf(".tab") > 0) "table/$filename"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else null
//         }

//     input: 
//     set file('blast.fmt7'), file('fasta.fasta'), val(id), val(source) from ch_igblast_filter
//     file imgtbase from ch_imgt_db_for_igblast_filter.mix(ch_imgt_db_for_igblast_filter_mix).collect()

//     output:
//     set source, id, file("${id}_parse-select.tab") into ch_for_merge
//     file "${id}_parse-select_sequences.fasta"
//     file "${id}_parse-select.tab"
//     file "${id}_command_log.txt" into igblast_log

//     when:
//     !params.downstream_only

//     script:
//     if (params.loci == 'ig'){
//         """
//         MakeDb.py igblast -i blast.fmt7 -s fasta.fasta -r \\
//         ${imgtbase}/${params.species}/vdj/imgt_${params.species}_IGHV.fasta \\
//         ${imgtbase}/${params.species}/vdj/imgt_${params.species}_IGHD.fasta \\
//         ${imgtbase}/${params.species}/vdj/imgt_${params.species}_IGHJ.fasta \\
//         --regions --scores > "${id}_command_log.txt"
//         ParseDb.py split -d blast_db-pass.tab -f FUNCTIONAL >> "${id}_command_log.txt"
//         ParseDb.py select -d blast_db-pass_FUNCTIONAL-T.tab -f V_CALL J_CALL -u "IGH" --regex --logic all --outname ${id} >> "${id}_command_log.txt"
//         ConvertDb.py fasta -d ${id}_parse-select.tab --if SEQUENCE_ID --sf SEQUENCE_IMGT --mf V_CALL DUPCOUNT >> "${id}_command_log.txt"
//         """
//     } else if (params.loci == 'tr') {
//         """
//         MakeDb.py igblast -i blast.fmt7 -s fasta.fasta -r \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRAV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRAJ.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRBV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRBD.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRBJ.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRDV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRDD.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRDJ.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRGV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRGJ.fasta" \\
//         --regions --scores > "${id}_command_log.txt"
//         ParseDb.py split -d blast_db-pass.tab -f FUNCTIONAL >> "${id}_command_log.txt"
//         ParseDb.py select -d blast_db-pass_FUNCTIONAL-T.tab -f V_CALL J_CALL -u "TR" --regex --logic all --outname ${id} >> "${id}_command_log.txt"
//         ConvertDb.py fasta -d ${id}_parse-select.tab --if SEQUENCE_ID --sf SEQUENCE_IMGT --mf V_CALL DUPCOUNT >> "${id}_command_log.txt"
//         """
//     }
// }

// //Merge tables belonging to the same patient
// process merge_tables{
//     tag "merge tables"
//     publishDir "${params.outdir}/genotyping/$source", mode: params.publish_dir_mode

//     input:
//     set source, id, file(tab) from ch_for_merge.groupTuple()

//     output:
//     set source, file("${source}.tab") into ch_for_shazam

//     when:
//     !params.downstream_only

//     script:
//     """
//     echo "${source}"
//     echo "${tab}"
//     echo "${tab[0]}"
//     echo "${tab.join('\n')}" > tab.list
    
//     head -n 1 ${tab[0]} > ${source}.tab
//     tail -n +2 ${tab} >> ${source}.tab
//     sed -i '/==>/d' ${source}.tab
//     """
// }


// if (params.loci == "ig"){
//     //Shazam! 
//     process shazam_ig{
//         tag "${id}"    
//         publishDir "${params.outdir}/genotyping/$id", mode: params.publish_dir_mode,
//             saveAs: {filename ->
//                 if (filename == "igh_genotyped.tab") "${id}_igh_genotyped.tab"
//                 else if (filename.indexOf("command_log.txt") > 0) "$filename"
//                 else if (filename == "threshold.txt" && !params.set_cluster_threshold) "${id}_threshold.txt"
//                 else if (filename == "genotype.pdf") "${id}_genotype.pdf"
//                 else if (filename == "Hamming_distance_threshold.pdf") "${id}_hamming_distance_threshold.pdf"
//                 else if (filename == "v_genotype.fasta") "${id}_v_genotype.fasta"
//                 else null
//             }

//         input:
//         set val(id), file(tab) from ch_for_shazam
//         file imgtbase from ch_imgt_db_for_shazam.mix(ch_imgt_db_for_shazam_mix).collect()

//         output:
//         set file("threshold.txt"), file("v_genotyped.tab"), val("$id") into ch_threshold_for_clone_definition_ig
//         file("v_genotype.fasta") into ch_fasta_for_clone_definition_ig
//         file "Hamming_distance_threshold.pdf" 
//         file "genotype.pdf"

//         when:
//         !params.downstream_only

//         script:
//         """
//         TIgGER-shazam.R $tab $params.loci $params.threshold_method ${imgtbase}/${params.species}/vdj/imgt_human_IGHV.fasta 
//         """

//     }

//     ch_threshold_for_clone_definition_tr = Channel.empty()
//     ch_fasta_for_clone_definition_tr = Channel.empty()
//     create_germlines_log_tr = Channel.empty()

// } else if (params.loci == "tr"){
//     //Shazam! 
//     process shazam_tr{
//         tag "${id}"    
//         publishDir "${params.outdir}/genotyping/$id", mode: params.publish_dir_mode,
//             saveAs: {filename ->
//                 if (filename == "v_tr_nogenotyped.tab") "${id}_igh_genotyped.tab"
//                 else if (filename.indexOf("command_log.txt") > 0) "$filename"
//                 else if (filename == "threshold.txt" && !params.set_cluster_threshold) "${id}_threshold.txt"
//                 else if (filename == "genotype.pdf") "${id}_genotype.pdf"
//                 else if (filename == "Hamming_distance_threshold.pdf") "${id}_hamming_distance_threshold.pdf"
//                 else if (filename == "v_genotype.fasta") "${id}_v_genotype.fasta"
//                 else null
//             }

//         input:
//         set val(id), file(tab) from ch_for_shazam
//         file imgtbase from ch_imgt_db_for_shazam.mix(ch_imgt_db_for_shazam_mix).collect()

//         output:
//         set file("threshold.txt"), file("v_tr_nogenotyped.tab"), val("$id") into ch_threshold_for_clone_definition_tr
//         file "Hamming_distance_threshold.pdf" 

//         when:
//         !params.downstream_only

//         script:
//         """
//         TIgGER-shazam.R $tab $params.loci $params.threshold_method \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRAV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRBV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRDV.fasta" \\
//         "${imgtbase}/${params.species}/vdj/imgt_${params.species}_TRGV.fasta"
//         """
//     }

//     ch_threshold_for_clone_definition_ig = Channel.empty()
//     ch_fasta_for_clone_definition_tr = Channel.from('no_file')
//     ch_fasta_for_clone_definition_ig = Channel.empty()
//     create_germlines_log_tr = Channel.empty()


// }

// //Assign clones
// process assign_clones{
//     tag "${id}" 
//     publishDir "${params.outdir}/clone_assignment/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf("table.tab") > 0) "${id}_log_table.tab"
//             else if (filename.indexOf("clone-pass.tab") > 0) "${id}_clone_pass.tab"
//             else if (filename.indexOf("command_log.txt") > 0) "$filename"
//             else if (filename == "threshold.txt" && !params.set_cluster_threshold) "${id}_threshold.txt"
//             else null
//         }

//     input:
//     set val(threshold), file(geno), val(id) from ch_threshold_for_clone_definition_ig.mix(ch_threshold_for_clone_definition_tr)
//     file(geno_fasta) from ch_fasta_for_clone_definition_ig.mix(ch_fasta_for_clone_definition_tr.collect())

//     output:
//     set file("${geno.baseName}_clone-pass.tab"), val("$id") into ch_for_germlines
//     file("$geno_fasta") into ch_fasta_for_germlines
//     file "${geno.baseName}_clone-pass.tab"
//     file "${geno.baseName}_table.tab"
//     file "${id}_command_log.txt" into assign_clones_log

//     when:
//     !params.downstream_only

//     script:
//     if (params.set_cluster_threshold) {
//         thr = params.cluster_threshold
//     } else {
//         thr = file(threshold).text
//         thr = thr.trim()
//     }
//     """
//     DefineClones.py -d $geno --act set --model ham --norm len --dist $thr --outname ${geno.baseName} --log ${geno.baseName}.log > "${id}_command_log.txt"
//     ParseLog.py -l "${geno.baseName}.log" -f ID VCALL JCALL JUNCLEN CLONED FILTERED CLONES
//     """
// }


// //Reconstruct germline sequences
// process germline_sequences{
//     tag "${id}"
//     publishDir "${params.outdir}/germlines/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf(".fasta") > 0) "fasta/$filename"
//             else if (filename.indexOf(".log") > 0) null
//             else if (filename.indexOf("table.tab") > 0) "$filename"
//             else if (filename.indexOf(".tab") > 0) "$filename"
//             else if (filename.indexOf("command_log.txt") >0) "$filename"
//             else null
//         }

//     input: 
//     set file(clones), val(id) from ch_for_germlines
//     file(geno_fasta) from ch_fasta_for_germlines
//     file imgtbase from ch_imgt_db_for_germline_sequences.mix(ch_imgt_db_for_germline_sequences_mix).collect()

//     output:
//     set val("$id"), file("${id}.tab") into ch_for_lineage_reconstruction
//     file "${id}.tab"
//     file "${id}_command_log.txt" into create_germlines_log_ig

//     when:
//     !params.downstream_only
//     params.loci == "ig"

//     script:
//     """
//     CreateGermlines.py -d ${clones} -g dmask --cloned -r $geno_fasta ${imgtbase}/human/vdj/imgt_human_IGHD.fasta ${imgtbase}/human/vdj/imgt_human_IGHJ.fasta --log ${clones.baseName}.log -o "${id}.tab" > "${id}_command_log.txt"
//     ParseLog.py -l "${clones.baseName}.log" -f ID V_CALL D_CALL J_CALL
//     """
// }



// //Lineage reconstruction
// process lineage_reconstruction{
//     tag "${id}"
//     publishDir "${params.outdir}/clonal_analysis/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf(".graphml") > 0) "$filename"
//             else if (filename.indexOf(".tsv") > 0) "$filename"
//             else if (filename.indexOf(".svg") > 0) "clone_tree_plots/$filename"
//             else null
//         }
    
//     input:
//     set val(id), file(clones) from ch_for_lineage_reconstruction

//     output:
//     set val("$id"), file("${id}.tab") into ch_for_clonal_analysis
//     file "lineage_reconstruction/*.tsv"
//     file "lineage_reconstruction/Clone_tree_plots/*.svg"
//     file "lineage_reconstruction/Graphml_trees/All_graphs_patient.graphml"

//     when:
//     !params.downstream_only
//     params.loci == "ig"

//     script:
//     """
//     which dnapars > dnapars_exec.txt
//     lineage_reconstruction.R
//     merge_graphs.sh
//     """
// }

// //Clonal analysis
// process clonal_analysis{
//     tag "${id}"
//     publishDir "${params.outdir}/clonal_analysis/$id", mode: params.publish_dir_mode,
//         saveAs: {filename ->
//             if (filename.indexOf(".tab") > 0) "$filename"
//             else if (filename.indexOf(".zip") > 0) "$filename"
//             else null
//         }
    
//     input:
//     set val(id), file(clones) from ch_for_clonal_analysis.mix(ch_tabs_for_clonal_analysis)

//     output:
//     file("${id}.tab") into ch_for_repertoire_comparison
//     file "clonal_analysis.zip"

//     when:
//     !params.skip_downstream

//     script:
//     """
//     clonal_analysis.R
//     zip -r clonal_analysis.zip clonal_analysis
//     """

// }

// //Repertoire comparison
// process repertoire_comparison{
//     tag "all" 
//     publishDir "${params.outdir}/repertoire_comparison", mode: params.publish_dir_mode,
//     saveAs: {filename ->
//             if (filename.indexOf(".tab") > 0) "table/$filename"
//             else if (filename.indexOf(".zip") > 0) "$filename"
//             else null
//         }


//     input:
//     file '*.tab' from ch_for_repertoire_comparison.collect()

//     output:
//     file '*.tab'
//     file "repertoire_comparison.zip"

//     when:
//     !params.skip_downstream

//     script:
//     """
//     repertoire_comparison.R
//     zip -r repertoire_comparison.zip repertoire_comparison
//     """
// }

// //Processing logs
// process processing_logs{
//     publishDir "${params.outdir}/parsing_logs", mode: params.publish_dir_mode

//     input:
//     file('filter_by_sequence_quality/*') from filter_by_sequence_quality_log.collect()
//     file('mask_primers/*') from mask_primers_log.collect()
//     file('pair_sequences/*') from pair_seq_log.collect()
//     file('cluster_sets/*') from cluster_sets_log.collect()
//     file('build_consensus/*') from build_consensus_log.collect()
//     file('repair_mates/*') from repair_log.collect()
//     file('assemble_pairs/*') from assemble_log.collect()
//     file('deduplicates/*') from dedup_log.collect()
//     file('filter_representative_2/*') from filter_seqs_log.collect()
//     file('igblast/*') from igblast_log.collect()
//     file('define_clones/*') from assign_clones_log.collect()
//     file('create_germlines/*') from create_germlines_log_ig.mix(create_germlines_log_tr).collect()
//     file('metadata.tsv') from ch_metadata_file_for_process_logs

//     output:
//     file "Table_sequences_process.tsv"

//     script:
//     '''
//     log_parsing.py
//     '''
// }

// //MultiQC
// process multiqc {
//     publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

//     input:
//     file (multiqc_config) from ch_multiqc_config
//     file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
//     // TODO nf-core: Add in log files from your new processes for MultiQC to find!
//     file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
//     file ('software_versions/*') from ch_software_versions_yaml.collect()
//     file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

//     output:
//     file "*multiqc_report.html" into ch_multiqc_report
//     file "*_data"
//     file "multiqc_plots"

//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
//     """
//     multiqc -f $rtitle $rfilename $custom_config_file .
//     """
// }
