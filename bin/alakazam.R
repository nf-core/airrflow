library(igraph)
#library(dplyr)
library(alakazam)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("Two arguments must be supplied (input file tab and input file fasta).\n", call.=FALSE)
} 

inputtable = args[1]
# "docker/QMKMK072AD/output_changeo_clonal/3463pbl_S8_L001_clone-pass_germ-pass.tab"

output_folder = dirname(inputtable)

db <- readChangeoDb(inputtable)

dbsubset <- subset(db, C_PRIMER %in% c("IgG", "IgA", "IgM"))

# Calculate and plot rank abundance curve
a <- estimateAbundance(dbsubset, group="C_PRIMER")
ggsave(paste(output_folder,"abundance.pdf",sep="/"), plotAbundanceCurve(a,silent=T))

# Generate Hill diversity curve
d <- rarefyDiversity(dbsubset, group="C_PRIMER")
ggsave(paste(output_folder,"diversity.pdf",sep="/"), plotDiversityCurve(d,silent=T))

# # Calculate CD3 amino acid properties
# p <- aminoAcidProperties(dbsubset, seq="JUNCTION", nt = T, trim=T, label="CDR3")
# 
# # V family usage by isotype and clone
# v <- countGenes(dbsubset, gene="V_CALL_GENOTYPED", groups="PRCONS", clone="CLONE", mode="family")
# 
# 
# # Example build tree for a single clone
# x <- makeChangeoClone(subset(db, CLONE == 37), text_fields="C_PRIMER")
# g <- buildPhylipLineage(x, dnapars="/Users/gisela/bin/phylip-3.697/exe/dnapars")
# 
# 
# #sample_clones <- readChangeoDb(inputfile)
# sample_clones <- db
# clones <- sample_clones %>%
#   group_by(CLONE) %>%
#   do(CHANGEO=makeChangeoClone(., text_fields="C_PRIMER",
#                               num_fields="DUPCOUNT"))
# 
# # Build lineages
# dnapars_exec <- "/Users/gisela/bin/phylip-3.697/exe/dnapars"
# graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
#                  dnapars_exec=dnapars_exec, rm_temp=TRUE)
# # Note, clones with only a single sequence will not be processed.
# # A warning will be generated and NULL will be returned by buildPhylipLineage
# # These entries may be removed for clarity
# graphs[sapply(graphs, is.null)] <- NULL
# graph <- graphs[[7]]
# plot(graph)
# # The graph has shared annotations for the clone
# data.frame(CLONE=graph$clone,
#            JUNCTION_LENGTH=graph$junc_len,
#            V_GENE=graph$v_gene,
#            J_GENE=graph$j_gene)
# 
# # The vertices have sequence specific annotations
# data.frame(SEQUENCE_ID=V(graph)$name, 
#            ISOTYPE=V(graph)$C_PRIMER,
#            DUPCOUNT=V(graph)$DUPCOUNT)
# 
# # Modifying layout
# V(graph)$color <- "steelblue"
# V(graph)$color[V(graph)$name == "Germline"] <- "black" 
# V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white" 
# V(graph)$label <- V(graph)$C_PRIMER
# E(graph)$label <- ""
# # Remove large default margins
# par(mar=c(0, 0, 0, 0) + 0.1)
# # Plot graph
# plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
#      vertex.label.color="black", vertex.size=20) # Add legend
# legend("topleft", c("Germline", "Inferred", "Sample"), fill=c("black", "white", "steelblue"), cex=0.75)
