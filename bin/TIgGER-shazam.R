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

IGHV_fasta = args[2]

output_folder = dirname(args[1])
  
db <- readChangeoDb(inputtable)
ighv <- readIgFasta(IGHV_fasta, strip_down_name = TRUE, force_caps = TRUE)

gt <- inferGenotype(db, germline_db = ighv, find_unmutated = FALSE)

gtseq <- genotypeFasta(gt, ighv)
writeFasta(gtseq, paste(output_folder,"v_genotype.fasta",sep="/"))

# Plot genotype
ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

# Modify allele calls and output TSV file
db_reassigned <- reassignAlleles(db, gtseq)
writeChangeoDb(db_reassigned, paste(output_folder,"igh_genotyped.tab",sep="/"))

################
#### shazam ####
################

dist_ham <- distToNearest(db_reassigned, vCallColumn="V_CALL_GENOTYPED", model="ham", 
                          normalize="len", nproc=1, first = FALSE)

# Find threshold using density method
output <- findThreshold(dist_ham$DIST_NEAREST, method="density")
threshold <- output@threshold


# Plot distance histogram, density estimate and optimum threshold
ggsave(paste(output_folder,"Hamming_distance_threshold.pdf",sep="/"), plot(output), device="pdf")

print(threshold)
write.table(threshold, file= paste(output_folder,"threshold.txt",sep="/"), quote=FALSE, sep="", row.names = FALSE, col.names = FALSE)
