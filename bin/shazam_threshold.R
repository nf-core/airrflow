#!/usr/bin/env Rscript
library(ggplot2)
library(stringi)
library(alakazam)
library(shazam)
library(dplyr)

# Set random seed for reproducibility
set.seed(12345)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
    stop("Two arguments must be supplied (input file tab and input file fasta).\n", call.=FALSE)
}

# Get input from args
inputtable = args[1]
threshold_method = args[2]


output_folder = dirname(inputtable)

db <- read.table(inputtable, header=TRUE, sep="\t")

# Add label for source species
sourceLabel <- gsub(pattern = "\\.tsv$", "", inputtable)

# Find the Hamming distance
dist_ham <- distToNearest(db,
                            vCallColumn="v_call",
                            jCallColumn="j_call",
                            sequenceColumn="junction",
                            model="ham",
                            normalize="len",
                            nproc=1,
                            first = FALSE)


num_dist <- length(unique(na.omit(dist_ham$dist_nearest)))

if (num_dist > 3) {
    # Find threshold using chosen method
    if (threshold_method == "density") {
        output <- findThreshold(dist_ham$dist_nearest, method="density")
        threshold <- output@threshold
    } else if (threshold_method == "gmm") {
        output <- findThreshold(dist_ham$dist_nearest, method="gmm")
        threshold <- output@threshold
    } else {
        stop("Threshold method is not available, please choose from: density, gmm")
    }

    # Plot distance histogram, density estimate and optimum threshold
    ggsave(paste(output_folder,paste0(sourceLabel, "_Hamming_distance_threshold.pdf"),sep="/"), plot(output), device="pdf")

} else {
    # Workaround for sources with too few nearest distance values to determine an effective threshold.
    # Set threshold to 0 and print a warning
    threshold <- 0.0
    warning(paste("Could not determine an effective Hamming distance threshold for source:", sourceLabel, ", which has", num_dist, "unique nearest distances. Threshold defaulting to 0.",  sep=" "))
    ggsave(paste(output_folder,paste0(sourceLabel, "_Hamming_distance_threshold.pdf"),sep="/"), plot(dist_ham$dist_nearest, dist_ham$duplicate_count), device="pdf")
}

write.table(threshold, file= paste(output_folder,paste0(sourceLabel, "_threshold.txt"),sep="/"), quote=FALSE, sep="", row.names = FALSE, col.names = FALSE)
