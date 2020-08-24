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

inputtable = args[1]

loci = args[2]

fastas = args[3:]

output_folder = dirname(inputtable)
  
db <- readChangeoDb(inputtable)

if (loci == "ig"){

  db_fasta <- readIgFasta(fastas, strip_down_name = TRUE, force_caps = TRUE)

  gt <- inferGenotype(db, find_unmutated = FALSE)

  gtseq <- genotypeFasta(gt, db_fasta)
  writeFasta(gtseq, paste(output_folder,"v_genotype.fasta",sep="/"))

  # Plot genotype
  ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

  # Modify allele calls and output TSV file
  db_reassigned <- reassignAlleles(db, gtseq)

  # Find the Hamming distance
  dist_ham <- distToNearest(db_reassigned, vCallColumn="V_CALL_GENOTYPED", model="ham", 
                          normalize="len", nproc=1, first = FALSE)
} else if (loci == "tr") {

  db_fasta_TRAV = fastas[1]
  db_fasta_TRBV = fastas[2]
  db_fasta_TRGV = fastas[3]
  db_fasta_TRGV = fastas[4]

  gt <- inferGenotype(db, find_unmutated = FALSE)

  gtseq <- genotypeFasta(gt, c(db_fasta_TRAV, db_fasta_TRBV, db_fasta_TRGV, db_fasta_TRGV))
  writeFasta(gtseq, paste(output_folder,"v_genotype.fasta",sep="/"))

  # Plot genotype
  ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

  # Modify allele calls and output TSV file
  db_reassigned <- reassignAlleles(db, gtseq)

  # Find the Hamming distance
  dist_ham <- distToNearest(db_reassigned, vCallColumn="V_CALL_GENOTYPED", model="ham", 
                          normalize="len", nproc=1, first = FALSE)
} else {
  stop("Loci specified is not available, please choose from: ig, tr.")
}

writeChangeoDb(db_reassigned, paste(output_folder,"v_genotyped.tab",sep="/"))

# Find threshold using density method
output <- findThreshold(dist_ham$DIST_NEAREST, method="density")
threshold <- output@threshold


# Plot distance histogram, density estimate and optimum threshold
ggsave(paste(output_folder,"Hamming_distance_threshold.pdf",sep="/"), plot(output), device="pdf")

print(threshold)
write.table(threshold, file= paste(output_folder,"threshold.txt",sep="/"), quote=FALSE, sep="", row.names = FALSE, col.names = FALSE)
