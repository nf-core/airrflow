#!/usr/bin/env Rscript
# 
# Chimeric sequences filter:
# sequences with >5 mismatches from the germline reference in any 10 base pair window are discarded
# Arguments:
#   --repertoire  Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
#   --collapseby Names of the metadata columns to be used to collapse
#              duplicated sequences. Default: NULL
#   -h  Display help.
# Example: ./collapseDuplicates.R --repertoire test-results/assign-genes/sc5p_v2_hs_PBMC_1k_b_airr_rearrangement_sequences_igblast_db-pass.tsv --ids input_id_1

# Libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("alakazam"))

  
# Define commmandline arguments
opt_list <- list(make_option(c("--repertoire"), dest="REPERTOIRE", default=NULL,
                             help="Input repertoire .tsv file"),
                 make_option(c("--collapseby"), dest="COLLAPSEBY", default=c("filename"),
                             help="Grouping fields to collapse duplicated sequences."),         
                 make_option(c("--single-cell"), dest="SINGLECELL", default=c("false"),
                             help="If true, cell_id will be added to the grouping fields to collapse duplicated sequences."),         
                 make_option(c("--ids"), dest="IDS", default=NULL,
                             help="Ids"),                                                            
                 make_option(c("--outname"), dest="OUTNAME", default=NULL,
                             help="Output name"))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))
opt

getwd()

# Check input file
if (!("REPERTOIRE" %in% names(opt))) {
    stop("You must provide an input repertoire file with --repertoire")
}

repertoires <- paste0(strsplit(opt$REPERTOIRE,".tsv,*")[[1]], ".tsv")
if (!is.null(opt$IDS)) {
    ids <- strsplit(opt$IDS,",")[[1]]
    if (length(ids) != length(repertoires)) {
        stop ("Number of repertoires and number of ids provided, don't match.")
    }
} else {
    stop("Error. Expecting repertoire ids. Provide --ids.")
}

# Read repertoire file
db <- bind_rows(lapply(repertoires, read_rearrangement))

num_fields <- c("consensus_count", "duplicate_count")
num_fields <- intersect(num_fields, colnames(db))

if (length(num_fields)<1) {
   num_fields <- NULL
}

collapseby <- strsplit(opt$COLLAPSEBY,",")[[1]]
singlecell <- strsplit(opt$SINGLECELL,",")[[1]]
if (singlecell) {
    message("Adding cell_id to the collapsing fields.")
    collapseby <- unique(c(collapseby,'cell_id'))
}
check <- checkColumns(db, collapseby)
if (check != TRUE) { stop(check) }

collapse_groups <- c("v_gene", 
                     "j_gene", 
                     "junction_length", 
                     "c_call", 
                     "productive",
                     collapseby)

db <- db %>%
   mutate(v_gene=getGene(v_call),
          j_gene=getGene(j_call)) %>%
   group_by(.dots=collapse_groups) %>%
   do(collapseDuplicates(.,
                         id = "sequence_id",
                         seq = "sequence_alignment",
                         text_fields = NULL,
                         num_fields = num_fields,
                         seq_fields = NULL,
                         add_count = TRUE,
                         ignore = c("N", "-", ".", "?"),
                         sep = ",",
                         dry = FALSE,
                         verbose = FALSE
   )) %>%
   ungroup() %>%
   select(-v_gene, -j_gene)


for (i in 1:length(repertoires)) {
    if (!is.null(opt$OUTNAME)) {   
        output_fn <- paste0(ids[i],"_",opt$OUTNAME,"_collapse-pass.tsv")
    } else {
       output_fn <- paste0(ids[i],"_collapse-pass.tsv")
    }
    write_airr(db %>% filter(id == ids[i]), file=output_fn)
}


