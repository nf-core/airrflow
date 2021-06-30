#!/usr/bin/env Rscript
# 
# Junction length multiple of 3 filter:
# Arguments:
#   --repertoire  Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
#   -h  Display help.
# Example: ./reveal_mod_3_junction.R --repertoire test-results/assign-genes/sc5p_v2_hs_PBMC_1k_b_airr_rearrangement_sequences_igblast_db-pass.tsv

# Libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("airr"))


# Define commmandline arguments
opt_list <- list(
    make_option(c("-r", "--repertoire"), dest="REPERTOIRE", default=NULL,help="Input repertoire .tsv file"),
    make_option(c("-o", "--outname"), dest="OUTNAME", default=NULL, help="Output name")
    )

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))
opt
getwd()

# Check input file
if (!("REPERTOIRE" %in% names(opt))) {
    stop("You must provide an input repertoire file with --repertoire")
}

# Read metadata file
db <- read_rearrangement(opt$REPERTOIRE)

# Filter and save
filter_pass <- db$junction_length %% 3 == 0 
table(filter_pass)

if (!is.null(opt$OUTNAME)) {
    output_fn <- paste0(opt$OUTNAME,"_junction-pass.tsv")
} else {
    output_fn <- sub(".tsv$", "_junction-pass.tsv", basename(opt$REPERTOIRE))
}

write_rearrangement(db[filter_pass,], file=output_fn)

