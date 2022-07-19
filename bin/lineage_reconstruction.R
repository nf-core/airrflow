#!/usr/bin/env Rscript
#
# Create lineage trees:
# Arguments:
#   --repertoire  Tabulated data in AIRR (TSV) format with clonal assignments and germline assignments.
#   --node-label Names of the metadata column to be used as node label on the tree plots
#   -h  Display help.
# Example: ./lineage_reconstruction.R --repertoire igblast_germ-pass.tsv --nodelabel population

# Libraries
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--node-label"), dest="node_text", default=c("filename"),
                help="Text to be used as node label. Provide 'none' if no label is desired.")
)

opt <- parse_args(OptionParser(option_list=opt_list))

theme_set(theme_bw(base_family = "ArialMT") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))

# Set output directories
patdir_lineage_trees <- "Clone_tree_plots"
dir.create(patdir_lineage_trees)
patdir_lineage_graphml <- "Graphml_trees"
dir.create(patdir_lineage_graphml)

# Read patient table
df_pat <- read.csv(opt$repertoire, sep="\t")

print(paste0("Node text request: ", opt$node_text))

avail_text = colnames(df_pat)

if (opt$node_text %in% avail_text) {
    print(paste0("Node string set to: ",opt$node_text))
} else if (opt$node_text == "none") {
    print("Node string set to: none.")
} else {
    print("Available fields: ")
    print(avail_text)
    print("or 'none'.")
    stop("Node string must be one of the above fields.")
}

# save clonal table
countclones <- countClones(df_pat, clone="clone_id", copy="duplicate_count")
write.table(countclones, paste("Clones_table_patient_", df_pat$subject_id[1],"_",df_pat$pcr_target_locus[1],".tsv", sep=""), quote=F, sep="\t", row.names = F)

# Get dnapars exec path
dnapars_exec_tab <- read.csv("dnapars_exec.txt", header=F)
dnapars_exec <- as.character(dnapars_exec_tab[1,1])

# Create clonal tree per clone
save_graph <- function(df_pat, clone_num){
    print(paste0("Started processing clone:",clone_num))
    sub_db_clone = dplyr::filter(df_pat, clone_id == clone_num) %>%
                            dplyr::mutate(across(everything(),as.character)) %>%
                            dplyr::mutate(across(c(junction_length,duplicate_count), as.numeric))

    # Make changeo clone and Build Phylip Lineage
    if ( opt$node_text == "none" ) {
        clone <- makeChangeoClone(sub_db_clone, text_fields = c("c_primer", "subject_id",
                                                        "sample_id", "clone_id"),
                                        num_fields = "duplicate_count")
        graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp = T, verbose = F)
        V(graph)$label <- ""
    } else {
        sub_db_clone <- sapply(sub_db_clone[,opt$node_text], as.character)
        clone <- makeChangeoClone(sub_db_clone, text_fields = append(c("c_primer", "subject_id",
                                                            "sample_id", "clone_id"), opt$node_text),
                                            num_fields = "duplicate_count")
        V(graph)$label <- V(graph)[,opt$node_text]
    }

    #Modify graph and plot attributes
    V(graph)$color <- "steelblue"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"

    # Remove large default margins
    par(mar=c(0, 0, 0, 0) + 0.1)
    vsize = V(graph)$duplicate_count
    vsize[is.na(vsize)] <- 1

    # Save graph in graphML format
    write_graph(graph, file=paste(patdir_lineage_graphml, "/Graph_", clone@data$subject_id[1], "_", clone@data$pcr_target_locus[1], "_clone_id_", clone_num, ".txt", sep=""), format = c("graphml"))

    # Plot tree
    pdf(paste(patdir_lineage_trees,"/Clone_tree_", clone@data$subject_id[1], "_clone_id_", clone_num, ".pdf", sep=""))
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

