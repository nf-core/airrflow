#!/usr/bin/env Rscript
# Written by Susanna Marquez and released under the MIT license (2021).

# Quality filter:
# - locus should match v_call chain
# - locus should match c_call (TODO)
# - sequence_alignment min length informative positions 200
# - max 10% N nucleotides
#
# Arguments:
#   --repertoire  Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
#   -h  Display help.
# Example: ./reveal_filter_quality.R --repertoire test-results/assign-genes/sc5p_v2_hs_PBMC_1k_b_airr_rearrangement_sequences_igblast_db-pass.tsv

# Libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("stringi"))

# Define commmandline arguments
opt_list <- list(
    make_option(c("-r", "--repertoire"), dest = "REPERTOIRE", default = NULL, help = "Input repertoire .tsv file"),
    make_option(c("-o", "--outname"), dest = "OUTPUT", default = NULL, help = "Output name")
)
# Parse arguments
opt <- parse_args(OptionParser(option_list = opt_list))


# Check input file
if (!("REPERTOIRE" %in% names(opt))) {
    stop("You must provide an input repertoire file with --repertoire")
}

# Read metadata file
db <- read_rearrangement(opt$REPERTOIRE)

# Remove rows that have NA values in all of v_call, d_call and j_call (still there when directly calling IgBlast)
filter_na <- is.na(db$v_call) & is.na(db$d_call) & is.na(db$j_call)
db <- db[!filter_na,]

# locus field and locus obtained from v_call should match
if (packageVersion("alakazam") < "1.0.3") {
    getLocus <- function(segment_call, first = TRUE, collapse = TRUE,
                        strip_d = TRUE, omit_nl = FALSE, sep = ",") {
        locus_regex <- "((IG[HLK]|TR[ABGD]))"
        r <- getSegment(segment_call, locus_regex,
            first = first, collapse = collapse,
            strip_d = strip_d, omit_nl = omit_nl, sep = sep
        )

        return(r)
    }
}

# Concordant locus
same_locus <- getLocus(db[["v_call"]]) == db[["locus"]]

# Max 10% N
n_count <- stri_count(db$sequence_alignment, regex = "Nn")
positions_count <- stri_count(db$sequence_alignment, regex = "[^-.]")
not_0 <- n_count > 0

if (any(not_0)) {
    n_count[not_0] <- n_count[not_0] / positions_count[not_0]
}
low_n <- n_count <= 0.10

# Min length 200 nt
long_seq <- stri_count(db$sequence_alignment, regex = "[^-.Nn]") >= 200

log <- data.frame(
    "same_locus" = same_locus,
    "low_n" = low_n,
    "long_seq" = long_seq, stringsAsFactors = F
)

summary <- log %>%
    group_by(same_locus, low_n, long_seq) %>%
    summarize(n = n(), .groups = "drop_last")

# Filter and save
filter_pass <- same_locus & low_n & long_seq

if (!is.null(opt$OUTPUT)) {
    output_fn <- paste0(opt$OUTPUT, "_quality-pass.tsv")
} else {
    output_fn <- sub(".tsv$", "_quality-pass.tsv", basename(opt$REPERTOIRE))
}
write_rearrangement(db[filter_pass, ], file = output_fn)

# cat("     TOTAL_GROUPS> ", n_groups,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)

write("START> FilterQuality", stdout())
write(paste0("FILE> ", basename(opt$REPERTOIRE)), stdout())
write(paste0("OUTPUT> ", basename(output_fn)), stdout())
write(paste0("PASS> ", sum(filter_pass)), stdout())
write(paste0("FAIL> ", sum(!filter_pass) + sum(filter_na)), stdout())
