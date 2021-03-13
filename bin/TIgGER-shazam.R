#!/usr/bin/env Rscript
library(ggplot2)
library(stringi)
library(alakazam)
library(tigger)
library(shazam)

# Set random seed for reproducibility
set.seed(12345)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
    stop("Two arguments must be supplied (input file tab and input file fasta).\n", call.=FALSE)
  } 

# Get input from args
inputtable = args[1]
loci = args[2]
threshold_method = args[3]
fastas = args[4:length(args)]

output_folder = dirname(inputtable)

#db <- readChangeoDb(inputtable)
db <- read.table(inputtable, header=TRUE, sep="\t")
print(colnames(db))

if (loci == "ig"){

  db_fasta <- readIgFasta(fastas, strip_down_name = TRUE)

  gt <- inferGenotype(db, v_call = "v_call", find_unmutated = T)

  gtseq <- genotypeFasta(gt, db_fasta)
  writeFasta(gtseq, paste(output_folder,"v_genotype.fasta",sep="/"))

  # Plot genotype
  ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

  # Modify allele calls and output TSV file
  db_reassigned <- reassignAlleles(db, gtseq)

  # Find the Hamming distance
  dist_ham <- distToNearest(db_reassigned, 
                            vCallColumn="v_call_genotyped",
                            jCallColumn="j_call",
                            sequenceColumn="junction",
                            model="ham", 
                            normalize="len", 
                            nproc=1, 
                            first = FALSE)
  writeChangeoDb(db_reassigned, paste(output_folder,"v_genotyped.tab",sep="/"))

} else if (loci == "tr") {

  db_fasta_TRAV = readIgFasta(fastas[1], strip_down_name = TRUE, force_caps = TRUE)
  db_fasta_TRBV = readIgFasta(fastas[2], strip_down_name = TRUE, force_caps = TRUE)
  db_fasta_TRDV = readIgFasta(fastas[3], strip_down_name = TRUE, force_caps = TRUE)
  db_fasta_TRGV = readIgFasta(fastas[4], strip_down_name = TRUE, force_caps = TRUE)

  gt <- inferGenotype(db, v_call = "v_call", find_unmutated = FALSE)

  gtseq <- genotypeFasta(gt, c(db_fasta_TRAV,db_fasta_TRBV,db_fasta_TRDV))
  writeFasta(gtseq, paste(output_folder,"TRxV_genotype.fasta",sep="/"))

  # Plot genotype
  ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

  # Modify allele calls and output TSV file
  db_reassigned <- reassignAlleles(db, gtseq)

  # Find the Hamming distance
  # TODO: check that this is fine (what about TRBV, etc.)
  dist_ham <- distToNearest(db_reassigned, 
                            vCallColumn="v_call",
                            jCallColumn="j_call",
                            sequenceColumn="junction",
                            model="ham", 
                            normalize="len", 
                            nproc=1, 
                            first = FALSE)
  
  writeChangeoDb(db, paste(output_folder,"v_tr_genotyped.tab",sep="/"))

} else {
  stop("Loci specified is not available, please choose from: ig, tr.")
}

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
ggsave(paste(output_folder,"Hamming_distance_threshold.pdf",sep="/"), plot(output), device="pdf")

print(threshold)
write.table(threshold, file= paste(output_folder,"threshold.txt",sep="/"), quote=FALSE, sep="", row.names = FALSE, col.names = FALSE)
