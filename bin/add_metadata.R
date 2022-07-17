#!/usr/bin/env Rscript
#
# Add metadata to repertoire:
# Arguments:
#   --repertoire    Tabulated data in AIRR (TSV) format with clonal assignments and germline assignments.
#   --samplesheet   Names of the metadata column to be used as node label on the tree plots
#   --outname       Filename for the output repertoire
#   -h  Display help.
# Example: ./lineage_reconstruction.R --repertoire igblast_germ-pass.tsv --nodelabel population

# Libraries
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dplyr))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--samplesheet"), default=NULL,
                help="Samplesheet providing extra metadata"),
    make_option(c("--outname"), default=NULL)
)

theme_set(theme_bw(base_family = "ArialMT") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

# Read repertoire
repertoire <- read.csv(opt$repertoire, sep="\t", header=TRUE, stringsAsFactors = F)
samplesheet <- read.cssv(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors = F)

parsed_fields <-
    c(
        "subject_id",
        "species",
        "pcr_target_locus",
        "filename_R1",
        "filename_R2",
        "filename_I1"
    )
metadata <- metadata[, !colnames(metadata) %in% parsed_fields]

# save repertoire table with metadata fields
anno_repertoire <- merge(x=repertoire, y=samplesheet, by.x = sample_id, by.y = sample_id, all.x=T)
write.table(anno_repertoire, opt$outname, quote=F, sep="\t", row.names = F)
