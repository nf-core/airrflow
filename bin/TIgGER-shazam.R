#!/usr/bin/env Rscript
library(ggplot2)
library(stringi)
library(alakazam)
library(tigger)
library(shazam)


set.seed(12345)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
    stop("Two arguments must be supplied (input file tab and input file fasta).\n", call.=FALSE)
  } 

inputtable = args[1]
#"docker/QMKMK072AD/output_changeo_blast/3463pbl_S8_L001_UMI_R1_R2_atleast-2_igblast_db-pass_FUNCTIONAL-T_parse-select.tab"

IGHV_fasta = args[2]
#"docker/imgt/human/vdj/imgt_human_IGHV.fasta"

output_folder = dirname(args[1])
  
db <- readChangeoDb(inputtable)
ighv <- readIgFasta(IGHV_fasta, strip_down_name = TRUE, force_caps = TRUE)

nv <- findNovelAlleles(db, germline_db = ighv, germline_min = 20)
gt <- inferGenotype(db, germline_db = ighv, find_unmutated = FALSE)

sel_nv <- selectNovel(nv)

if (nrow(sel_nv) > 0){
  ggsave(paste(output_folder,"novel_alleles.pdf",sep="/"), plotNovel(db, novel_row = sel_nv))
}

# Save genotype
gtseq <- genotypeFasta(gt, ighv, nv)
writeFasta(gtseq, paste(output_folder,"v_genotype.fasta",sep="/"))

# Plot genotype
ggsave(paste(output_folder,"genotype.pdf",sep="/"), plotGenotype(gt, silent=T))

# Modify allele calls and output TSV file
db_reassigned <- reassignAlleles(db, gtseq)
writeChangeoDb(db_reassigned, paste(output_folder,"igh_genotyped.tab",sep="/"))

#### shazam ####

dist_ham <- distToNearest(db_reassigned, vCallColumn="V_CALL_GENOTYPED", model="ham", 
                          normalize="len", nproc=1, first = FALSE)

# library(ggplot2)
# p1 <- ggplot(subset(dist_ham, !is.na(DIST_NEAREST)),
#              aes(x=DIST_NEAREST)) + 
#   theme_bw() + 
#   xlab("Hamming distance") + 
#   ylab("Count") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)

# Find threshold using density method
output <- findThreshold(dist_ham$DIST_NEAREST, method="density")
threshold <- output@threshold


# Plot distance histogram, density estimate and optimum threshold
ggsave(paste(output_folder,"Hamming_distance_threshold.pdf",sep="/"), plot(output), device="pdf")

print(threshold)
write.table(threshold, file= paste(output_folder,"threshold.txt",sep="/"), quote=FALSE, sep="", row.names = FALSE, col.names = FALSE)
