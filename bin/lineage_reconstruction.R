#!/usr/bin/env Rscript

library(alakazam)
library(ggplot2)
library(data.table)
library(dplyr)
library(tigger)
library(shazam)
library(igraph)
library(gplots)
library(circlize)
library(UpSetR)
library(gtools)

theme_set(theme_bw(base_family = "ArialMT") + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

datadir <- "."
patdir_lineage <- "lineage_reconstruction"
dir.create(patdir_lineage)
patdir_lineage_trees <- paste(patdir_lineage, "Clone_tree_plots", sep = "/")
dir.create(patdir_lineage_trees)
patdir_lineage_graphml <- paste(patdir_lineage, "Graphml_trees", sep = "/")
dir.create(patdir_lineage_graphml)

# Read patient table
fname <- system(paste0("find '",datadir,"' -name '*.tab'"), intern=T)
df_pat <- read.csv(fname, sep="\t")
df_pat$SAMPLE <- as.factor(paste(df_pat$TREATMENT, df_pat$EXTRACT_TIME, df_pat$SOURCE, sep="_"))
df_pat$SAMPLE_POP <- as.factor(paste(df_pat$TREATMENT, df_pat$EXTRACT_TIME, df_pat$SOURCE, df_pat$POPULATION, sep="_"))

##################
# Clonal lineages
#################

# save clonal table
countclones <- countClones(df_pat,clone="CLONE", copy="DUPCOUNT")
write.table(countclones, paste(patdir_lineage, "/", "Clones_table_patient_", df_pat$SOURCE[1],".tsv", sep=""), quote=F, sep="\t", row.names = F)

# Restrict clonal tree size
clones <- subset(countclones, SEQ_COUNT > 2 & SEQ_COUNT < 1000)

dnapars_exec_tab <- read.csv("dnapars_exec.txt", header=F)
dnapars_exec <- as.character(dnapars_exec_tab[1,1])

save_graph <- function(df_pat, clone_id){
    print(paste0("Started processing clone:",clone_id))
    sub_db_clone <- subset(df_pat, CLONE == clone_id)
    sub_db_clone$CLONE <- sapply(sub_db_clone$CLONE, as.character)
    sub_db_clone$SAMPLE <- sapply(sub_db_clone$SAMPLE, as.character)
    sub_db_clone$SAMPLE_POP <- sapply(sub_db_clone$SAMPLE_POP, as.character)
    sub_db_clone$C_PRIMER <- sapply(sub_db_clone$C_PRIMER, as.character)
    sub_db_clone$TREATMENT <- sapply(sub_db_clone$TREATMENT, as.character)
    sub_db_clone$POPULATION <- sapply(sub_db_clone$POPULATION, as.character)
    sub_db_clone$SOURCE <- sapply(sub_db_clone$SOURCE, as.character)
    sub_db_clone$EXTRACT_TIME <- sapply(sub_db_clone$EXTRACT_TIME, as.character)
    sub_db_clone$C_PRIMER <- sapply(sub_db_clone$C_PRIMER, as.character)
    
    clone <- makeChangeoClone(sub_db_clone, text_fields = c("C_PRIMER", "TREATMENT", "POPULATION", "SOURCE", "EXTRACT_TIME", "SAMPLE", "SAMPLE_POP", "CLONE"), num_fields = "DUPCOUNT")


    graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp = T)
    
    #Modify graph and plot attributes
    V(graph)$color <- "steelblue"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
    V(graph)$label <- V(graph)$POPULATION

    # Remove large default margins
    par(mar=c(0, 0, 0, 0) + 0.1)
    vsize = V(graph)$DUPCOUNT
    vsize[is.na(vsize)] <- 1

    # Save graph in graphML format
    write_graph(graph, file=paste(patdir_lineage_graphml, "/Graph_", clone@data$SOURCE[1],  "_clone_id_", clone_id, ".txt", sep=""), format = c("graphml"))

    # Plot
    svg(filename = paste(patdir_lineage_trees,"/Clone_tree_", clone@data$SOURCE[1], "_clone_id_", clone_id, ".svg", sep=""))
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
        vertex.label.color="black", vertex.size=(vsize/20 + 6))
    legend("topleft", c("Germline", "Inferred", "Sample"), 
        fill=c("black", "white", "steelblue"), cex=0.75)
    dev.off()
    
}

for (clone_id in clones$CLONE){
    tryCatch(withCallingHandlers(save_graph(df_pat, clone_id), 
                   error=function(e) {print(paste0("Skipping clone due to problem:", clone_id))},
                   warning=function(w) {print(paste0("Warning for clone:", clone_id))
                               invokeRestart("muffleWarning")}), 
           error = function(e) { print(paste0("Processed clone:", clone_id)) })
}

