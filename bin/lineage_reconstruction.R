#!/usr/bin/env Rscript

library(alakazam)
library(igraph)
library(dplyr)

theme_set(theme_bw(base_family = "ArialMT") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
    stop("Input file argument must be supplied.\n", call.=FALSE)
}

# Get input table from args
inputtable = args[1]
node_text = args[2]

print(paste0("Node text: ", node_text))

avail_text = c("c_primer", "treatment", "population", "source",
                "extract_time", "sample", "sample_pop", "clone_id", "seq_id", "none")

if (node_text %in% avail_text) {
    print(paste0("Node string set to: ",node_text))
} else {
    stop(paste0("Node string must be one of: ", avail_text))
}

# Set output directories
patdir_lineage_trees <- "Clone_tree_plots"
dir.create(patdir_lineage_trees)
patdir_lineage_graphml <- "Graphml_trees"
dir.create(patdir_lineage_graphml)

# Read patient table
df_pat <- read.csv(inputtable, sep="\t")
df_pat$sample <- as.factor(paste(df_pat$treatment, df_pat$extract_time, df_pat$source, sep="_"))
df_pat$sample_pop <- as.factor(paste(df_pat$treatment, df_pat$extract_time, df_pat$source, df_pat$population, sep="_"))

# save clonal table
countclones <- countClones(df_pat, clone="clone_id", copy="duplicate_count")
write.table(countclones, paste("Clones_table_patient_", df_pat$source[1],".tsv", sep=""), quote=F, sep="\t", row.names = F)

# Restrict clonal tree size
clones <- filter(countclones, seq_count > 2 & seq_count < 1000)
write.table(clones, paste("Clones_table_patient_filtered_between_3_and_1000_", df_pat$source[1],".tsv", sep=""), quote=F, sep="\t", row.names = F)

# Get dnapars exec path
dnapars_exec_tab <- read.csv("dnapars_exec.txt", header=F)
dnapars_exec <- as.character(dnapars_exec_tab[1,1])

# Create clonal tree per clone
save_graph <- function(df_pat, clone_num){
    print(paste0("Started processing clone:",clone_num))
    sub_db_clone <- subset(df_pat, clone_id == clone_num)
    sub_db_clone$clone_id <- sapply(sub_db_clone$clone_id, as.character)
    sub_db_clone$sample_id <- sapply(sub_db_clone$sample_id, as.character)
    sub_db_clone$c_primer <- sapply(sub_db_clone$c_primer, as.character)
    sub_db_clone$group_name <- sapply(sub_db_clone$group_name, as.character)
    sub_db_clone$subject_id <- sapply(sub_db_clone$subject_id, as.character)


    # Make changeo clone
    clone <- makeChangeoClone(sub_db_clone, text_fields = c("c_primer", "subject_id",
                                                            "sample_id", "clone_id", "group_name"),
                                            num_fields = "duplicate_count")

    # Build Phylip lineage
    graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp = T, verbose = F)

    #Modify graph and plot attributes
    V(graph)$color <- "steelblue"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"

    # Set label on the nodes
    if ( node_text == "c_primer" ) {
        V(graph)$label <- V(graph)$c_primer
    } else if ( node_text == "sample_id" ) {
        V(graph)$label <- V(graph)$sample_id
    } else if ( node_text == "clone_id" ) {
        V(graph)$label <- V(graph)$clone_id
    } else if ( node_text == "seq_id" ){
        V(graph)$label <- V(graph)$name
    } else if ( node_text == "none" ) {
        V(graph)$label <- ""
    }


    # Remove large default margins
    par(mar=c(0, 0, 0, 0) + 0.1)
    vsize = V(graph)$duplicate_count
    vsize[is.na(vsize)] <- 1

    # Save graph in graphML format
    write_graph(graph, file=paste(patdir_lineage_graphml, "/Graph_", clone@data$source[1],  "_clone_id_", clone_num, ".txt", sep=""), format = c("graphml"))

    # Plot tree
    pdf(paste(patdir_lineage_trees,"/Clone_tree_", clone@data$source[1], "_clone_id_", clone_num, ".pdf", sep=""))
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
        vertex.label.color="black", vertex.size=(vsize/20 + 6))
    legend("topleft", c("Germline", "Inferred", "Sample"),
        fill=c("black", "white", "steelblue"), cex=0.75)
    dev.off()

}

for (clone_num in clones$clone_id){
    tryCatch(withCallingHandlers(save_graph(df_pat, clone_num),
                    error=function(e) {print(paste0("Skipping clone due to problem:", clone_num))
                                        print("Here is the original error message:")
                                        print(e)},
                    warning=function(w) {print(paste0("Warning for clone:", clone_num))
                                        invokeRestart("muffleWarning")}),
                    error = function(e) { print(paste0("Processed clone:", clone_num)) })
}

