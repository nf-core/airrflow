#!/usr/bin/env Rscript
#
# Add metadata to repertoire:
# Arguments:
#   --repertoire    Tabulated data in AIRR (TSV) format with clonal assignments and germline assignments.
#   --samplesheet   Names of the metadata column to be used as node label on the tree plots
#   --outname       Filename for the output repertoire
#   -h  Display help.
# Example: ./add_metadata.R --repertoire igblast_germ-pass.tsv --samplesheet samplesheet.tsv --outname my-repertoire

# Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--samplesheet"), default=NULL,
                help="Samplesheet providing extra metadata"),
    make_option(c("--outname"), default=NULL)
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Read repertoire
repertoire <- read.csv(opt$repertoire, sep="\t", header=TRUE, stringsAsFactors = F)
samplesheet <- read.csv(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors = F)

parsed_fields <-
    c(
        "subject_id",
        "species",
        "pcr_target_locus",
        "filename_R1",
        "filename_R2",
        "filename_I1"
    )

samplesheet_colnames <- colnames(samplesheet)

# merge tables only in case the samplesheet contains more columns than the required ones
print( samplesheet_colnames[!(samplesheet_colnames %in% parsed_fields)])

if (length(samplesheet_colnames[!(samplesheet_colnames %in% parsed_fields)]) > 1 ) {
    print("None in parsed fields")
    samplesheet <- samplesheet[, !colnames(samplesheet) %in% parsed_fields]
    anno_repertoire <- base::merge(x=repertoire, y=samplesheet, by.x = "sample_id", by.y = "sample_id", all.x=T)
} else {
    print("Some in parsed fields")
    anno_repertoire <- repertoire
}

# save repertoire table with metadata fields
write.table(anno_repertoire, opt$outname, quote=F, sep="\t", row.names = F)
