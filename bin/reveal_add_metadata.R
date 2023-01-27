#!/usr/bin/env Rscript
#
# Merge db and input metadata by
#
# Arguments:
#   --repertoire Repertoire file
#   --input   Validated input metadata (contains input_id field)
#   --input_id input_id value in the `--input` file that will be used to
#              select the row that contains the metadata for the `--repertoire`
#              file
#   -h  Display help.
# Example: ./merge_db_input.R --input ../../test-datasets/metadata.tsv

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("dplyr"))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"),
        dest = "REPERTOIRE", default = NULL,
        help = "Repertoire tabulated file."
    ),
    make_option(c("--metadata"),
        dest = "METADATA", default = NULL,
        help = "Repertoire's validated metadata file."
    ),
    make_option(c("--input_id"),
        dest = "INPUTID", default = NULL,
        help = "input_id value in the `--metadata` file that
        will be used to select the row that contains the
        metadata for the `--repertoire` file"
    ),
    make_option(c("--outname"),
        dest = "OUTNAME",
        default = NULL, help = "Output name"
    )
)
# Parse arguments
opt <- parse_args(OptionParser(option_list = opt_list))
# opt
# getwd()

# Check repertoire file
if (!("REPERTOIRE" %in% names(opt))) {
    stop("You must provide a repertoire file with --repertoire.")
}

# Check input file
if (!("METADATA" %in% names(opt))) {
    stop("You must provide an input file with --metadata")
}

# Check input_id
if (!("INPUTID" %in% names(opt))) {
    stop("You must provide an input_id file with --input_id")
}

# Read metadata file
metadata <- read.csv(opt$METADATA, sep = "\t", header = TRUE, stringsAsFactors = F)

metadata <- metadata %>%
    filter(sample_id == opt$INPUTID)

if (nrow(metadata) != 1) {
    stop("Expecting nrow(metadata) == 1; nrow(metadata) == ", nrow(metadata), " found")
}

internal_fields <-
    c(
        "valid_filename",
        "valid_species",
        "valid_collapseby",
        "collapseby_group",
        "collapseby_size",
        "valid_cloneby",
        #        "cloneby_group",
        "cloneby_size",
        "id",
        "filetype",
        "valid_single_cell",
        "valid_pcr_target_locus",
        "filename_R1",
        "filename_R2",
        "filename_I1"
    )
metadata <- metadata[, !colnames(metadata) %in% internal_fields]

db <- read_rearrangement(opt$REPERTOIRE)

db <- cbind(db, metadata)

if (!is.null(opt$OUTNAME)) {
    output_fn <- paste0(opt$OUTNAME, "_meta-pass.tsv")
} else {
    output_fn <- sub(".tsv$", "_meta-pass.tsv", basename(opt$REPERTOIRE))
}

write_rearrangement(db, file = output_fn)


write("START> AddMetadata", stdout())
write(paste0("FILE> ", basename(opt$REPERTOIRE)), stdout())
write(paste0("OUTPUT> ", basename(output_fn)), stdout())
write(paste0("PASS> ", nrow(db)), stdout())
write(paste0("FAIL> ", 0), stdout())
