#!/usr/bin/env Rscript
#
# Chimeric sequences filter:
# sequences with >5 mismatches from the germline reference in any 10 base pair window are discarded
# Arguments:
#   --repertoire  Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
#   --outname     Output name
#   -h  Display help.
# Example: ./reveal_chimeric.R --repertoire test-results/assign-genes/db-pass.tsv

# Libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("shazam"))

mutThresh <- 6
windowSize <- 10

# Define commmandline arguments
opt_list <- list(make_option(c("--repertoire"), dest="REPERTOIRE", default=NULL,
            help="Input repertoire .tsv file"),
            make_option(c("--outname"), dest="OUTNAME", default=NULL,
            help="Output name"))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))
#opt
#getwd()

# Check input file
if (!("REPERTOIRE" %in% names(opt))) {
    stop("You must provide an input repertoire file with --repertoire")
}

# Read metadata file
db <- read_airr(opt$REPERTOIRE)

is_chimeric <- slideWindowDb(
    db,
    sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    mutThresh,
    windowSize
)

if (!is.null(opt$OUTNAME)) {
    output_fn <- paste0(opt$OUTNAME,"_chimera-pass.tsv")
} else {
    output_fn <- sub(".tsv$", "_chimera-pass.tsv", basename(opt$REPERTOIRE))
}

write_airr(db[!is_chimeric,] %>% select(-germline_alignment_d_mask), file=output_fn)

db$is_chimeric <- is_chimeric
write_airr(db %>% select(sequence_id, is_chimeric, sequence_alignment, germline_alignment_d_mask), file=sub(".tsv",".log.txt",output_fn))

write("START> RemoveChimeric", stdout())
write(paste0("FILE> ",basename(opt$REPERTOIRE)), stdout())
write(paste0("OUTPUT> ",basename(output_fn)), stdout())
write(paste0("PASS> ",sum(!is_chimeric)), stdout())
write(paste0("FAIL> ",sum(is_chimeric)), stdout())
