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
outdir <- "clonal_analysis"
dir.create(outdir)
patdir_overlap <- paste(outdir,"Clone_overlap",sep="/")
dir.create(patdir_overlap)
patdir_number <- paste(outdir,"Clone_numbers",sep="/")
dir.create(patdir_number)
patdir_lineage <- paste(outdir,"Clone_lineage",sep="/")
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


############################
# Clone chordplots
############################

# Chordplot comparison time points per patient
df_pat$TIME_POP <- as.factor(paste(df_pat$EXTRACT_TIME, df_pat$POPULATION, sep="_"))
df_pop_time <- split(df_pat, df_pat$TIME_POP)

count_clones <- countClones(df_pat)

## Calculating overlaps between time points
combin <- expand.grid(unique(df_pat$EXTRACT_TIME), unique(df_pat$POPULATION))
combin$names <- apply(combin[,c("Var1", "Var2")], 1, paste, collapse = "_")

baselines <- subset(combin, combin$Var1 == "baseline")
other <- subset(combin, combin$Var1 != "baseline")
clonedf <- expand.grid(baselines$names, other$names)
colnames(clonedf) <- c("from","to")
seqdf <- clonedf

lenintersects = numeric(0)
seqsintersects = numeric(0)
for (j in c(1:nrow(clonedf))){

    inter <- intersect(df_pop_time[[which(grepl(clonedf[j,1], names(df_pop_time)))]]$CLONE, 
                        df_pop_time[[which(grepl(clonedf[j,2], names(df_pop_time)))]]$CLONE)

    clones_subset <- count_clones[which(count_clones$CLONE %in% as.character(inter)),]

    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$SEQ_COUNT))
}

clonedf$value <- lenintersects
seqdf$value <- seqsintersects

# Saving both tables
write.table(clonedf, file = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(seqdf, file = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)

grid.col = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00")

# Clone overlap plot
svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".svg", sep=""))
chordDiagram(clonedf, grid.col = grid.col, self.link = 1,
            transparency = 0.3,
            annotationTrack="grid",
            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
            adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE OVERLAP", df_pop_time[[1]]$TREATMENT[1], df_pop_time[[1]]$SOURCE[1]), cex = 0.8)
circos.clear()
dev.off()

png(filename = paste(patdir_overlap, "/Clone_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
chordDiagram(clonedf, grid.col = grid.col, self.link = 1,
            transparency = 0.3,
            annotationTrack="grid",
            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
            adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE OVERLAP", df_pop_time[[1]]$TREATMENT[1], df_pop_time[[1]]$SOURCE[1]), cex = 0.8)
circos.clear()
dev.off()

# Sequences overlap plot
svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".svg", sep=""))
chordDiagram(seqdf, grid.col = grid.col, self.link = 1,
            transparency = 0.3,
            annotationTrack="grid",
            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
            adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$TREATMENT[1], df_pop_time[[1]]$SOURCE[1]), cex = 0.8)
circos.clear()
dev.off()

png(filename = paste(patdir_overlap, "/Clone_seqN_overlap_comparison_time_points_", df_pop_time[[1]]$TREATMENT[1], "_", df_pop_time[[1]]$SOURCE[1], ".png", sep=""), width=15, height=15, units = "cm", res = 300)
chordDiagram(seqdf, grid.col = grid.col, self.link = 1,
            transparency = 0.3,
            annotationTrack="grid",
            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
            adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE SEQ NUM OVERLAP", df_pop_time[[1]]$TREATMENT[1], df_pop_time[[1]]$SOURCE[1]), cex = 0.8)
circos.clear()
dev.off()



df_TP <- split(df_pat, df_pat$EXTRACT_TIME)

# Plots per patient and time point - overlap populations
#-------------------------------------------------------

for (n in c(1:length(df_TP))) {
df_pop <- split(df_TP[[n]], df_TP[[n]]$POPULATION)
vennplot <- venn(list(unique(df_pop[[1]]$CLONE), unique(df_pop[[2]]$CLONE), unique(df_pop[[3]]$CLONE), unique(df_pop[[4]]$CLONE)), names = names(df_pop))

listInput <- list(df_pop[[1]]$CLONE, df_pop[[2]]$CLONE, df_pop[[3]]$CLONE, df_pop[[4]]$CLONE)
names(listInput) <- names(df_pop)
combin <- data.frame(from=combinations(4,2,names(df_pop),repeats.allowed=F)[,1], to=combinations(4,2,names(df_pop),repeats.allowed=F)[,2])

listInput <- list(df_pop[[1]]$CLONE, df_pop[[2]]$CLONE, df_pop[[3]]$CLONE, df_pop[[4]]$CLONE)
names(listInput) <- names(df_pop)

svg(filename = paste(patdir_overlap,"/Set_plot_", df_pop[[1]]$TREATMENT[1], "_",df_pop[[1]]$EXTRACT_TIME[1], "_",df_pop[[1]]$SOURCE[1], ".svg", sep=""))
upset(fromList(listInput), group.by = "sets", order.by="freq", point.size = 3.5, line.size=2, mainbar.y.label = "Clone intersections", sets.x.label = "Clones per population")
dev.off()

png(filename = paste(patdir_overlap,"/Set_plot_", df_pop[[1]]$TREATMENT[1], "_",df_pop[[1]]$EXTRACT_TIME[1], "_",df_pop[[1]]$SOURCE[1], ".png", sep=""), res = 600, width = 15, height=10, units = "cm")
upset(fromList(listInput), order.by="freq", group.by = "sets", point.size = 3.5, line.size=2, mainbar.y.label = "Clone intersections", sets.x.label = "Clones per population")
dev.off()

clonedf <- combin
seqdf <- combin

lenintersects = numeric(0)
seqsintersects = numeric(0)
for (j in c(1:nrow(clonedf))){
    inter <- intersect(df_pop_time[[which(grepl(paste0("^",clonedf[j,1]), names(df_pop)))]]$CLONE, 
                        df_pop_time[[which(grepl(paste0("^",clonedf[j,2]), names(df_pop)))]]$CLONE)
    
    clones_subset <- count_clones[which(count_clones$CLONE %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$SEQ_COUNT))
}

clonedf$value <- lenintersects
seqdf$value <- seqsintersects


self_comb <- data.frame(from = names(df_pop), to = names(df_pop))
self_clonedf <- self_comb
self_seqdf<- self_comb

lenintersects <- numeric(0)
seqsintersects <- numeric(0)
for (pop in self_comb$from){
    inter <- attributes(vennplot)[["intersections"]][[pop]]
    clones_subset <- count_clones[which(count_clones$CLONE %in% as.character(inter)),]
    
    lenintersects <- c(lenintersects, length(inter))
    seqsintersects <- c(seqsintersects, sum(clones_subset$SEQ_COUNT))
}
self_clonedf$value <- lenintersects
self_seqdf$value <- seqsintersects

clonedf <- rbind(clonedf, self_clonedf)
seqdf <- rbind(seqdf, self_seqdf)

write.table(clonedf, file = paste(patdir_overlap,"/Clone_overlap_comparison_population_", df_pop[[1]]$TREATMENT[1], "_", df_pop[[1]]$EXTRACT_TIME[1], "_", df_pop[[1]]$SOURCE[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
write.table(seqdf, file = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_population_", df_pop[[1]]$TREATMENT[1], "_", df_pop[[1]]$EXTRACT_TIME[1], "_", df_pop[[1]]$SOURCE[1], ".tsv", sep=""), sep = "\t", quote = F, row.names = F)

grid.col = c("#225ea8","#41b6c4","#a1dab4","#ffffcc")


# Plots clone overlap
svg(filename = paste(patdir_overlap,"/Clone_overlap_comparison_population_", df_pop[[1]]$TREATMENT[1], "_", df_pop[[1]]$EXTRACT_TIME[1], "_", df_pop[[1]]$SOURCE[1], ".svg", sep=""))
chordDiagram(clonedf, grid.col = grid.col, self.link = 1,
                transparency = 0.3,
                annotationTrack="grid",
                preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(clonedf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE OVERLAP", df_pop[[1]]$TREATMENT[1], df_pop[[1]]$SOURCE[1], df_pop[[1]]$EXTRACT_TIME[1]), cex = 0.8)
circos.clear()
dev.off()

# Plots clone sequence numbers overlap
svg(filename = paste(patdir_overlap,"/Clone_seqN_overlap_comparison_population_", df_pop[[1]]$TREATMENT[1], "_", df_pop[[1]]$EXTRACT_TIME[1], "_", df_pop[[1]]$SOURCE[1], ".svg", sep=""))
chordDiagram(seqdf, grid.col = grid.col, self.link = 1,
            transparency = 0.3,
            annotationTrack="grid",
            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(seqdf))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
                adj = c(0, 0.5))
}, bg.border = NA)
title(paste("CLONE SEQ NUM OVERLAP", df_pop[[1]]$TREATMENT[1], df_pop[[1]]$SOURCE[1], df_pop[[1]]$EXTRACT_TIME[1]), cex = 0.8)
circos.clear()
dev.off()
}

#########################
## Clone numbers
#########################

## Clone numbers per sample
clones <- countClones(df_pat, groups=c("SAMPLE"))
clones$TREATMENT <- sapply(clones$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[1])
clones$TIME_POINT <- sapply(clones$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[3])
clones$PATIENT <- sapply(clones$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[2])


n_fun <- function(x){
return(data.frame(y = 2500, label = paste0("N = ",length(x))))
}

pclones1 <- ggplot(data = clones, aes(x=TIME_POINT, y=SEQ_COUNT)) +
geom_jitter() +
stat_summary(fun.y=mean, geom="point", color="red") +
stat_summary(fun.data=n_fun, geom="text") +
xlab("") +
ylab("N sequences / clone") +
facet_grid(cols=vars(PATIENT), drop=T)
ggsave(plot = pclones1, filename = paste(patdir_number,"/Clones_sequences_patient", "_", clones$PATIENT[1], "_", clones$TREATMENT[1],".svg",sep = ""), device="svg", width = 10, height = 7, units="cm")
ggsave(plot = pclones1, filename = paste(patdir_number,"/Clones_sequences_patient_", clones$PATIENT[1], "_", clones$TREATMENT[1],".png",sep = ""), device="png", width = 10, height = 7, units="cm")


clone_N <- count(clones, vars="SAMPLE")
clone_N$vars <- NULL

clone_mean <- tapply(clones$SEQ_COUNT, clones$SAMPLE, mean)
clone_mean <- data.frame(as.list(clone_mean))
clone_mean <- t(clone_mean)
clone_mean <- as.data.frame(clone_mean)
colnames(clone_mean) <- c("mean")
clone_mean$sample <- row.names(clone_mean)
row.names(clone_mean) <- NULL

clone_median <- tapply(clones$SEQ_COUNT, clones$SAMPLE, median)
clone_median <- data.frame(as.list(clone_median))
clone_median <- t(clone_median)
clone_median <- as.data.frame(clone_median)
colnames(clone_median) <- c("median")

clone_stats <- merge(x=clone_mean, y=clone_median, by.x = "sample", by.y = "row.names", all.x=T)
clone_number_stats <- merge(x = clone_N, y=clone_stats, by.x = "SAMPLE", by.y = "sample")
colnames(clone_number_stats) <- c("Sample", "Number of clones", "Mean number of sequences per clone", "Median number of sequences per clone")
write.table(clone_number_stats, file =paste(patdir_number,"/Clones_number_stats_patient_", clones$PATIENT[1], "_", clones$TREATMENT[1], ".tsv", sep=""), sep = "\t", row.names = F, quote = F)

# Clone numbers per population
clones <- countClones(df_pat, groups=c("SAMPLE_POP"))
clones$TREATMENT <- sapply(clones$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[1])
clones$TIME_POINT <- sapply(clones$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[2])
clones$TIME_POINT <- factor(clones$TIME_POINT, levels=c("baseline","6months"))
clones$PATIENT <- sapply(clones$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[3])
clones$POPULATION <- sapply(clones$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[4])


n_fun <- function(x){
return(data.frame(y = 2500, label = paste0("N = ",length(x))))
}

pclones2 <- ggplot(data = clones, aes(x=TIME_POINT, y=SEQ_COUNT)) +
geom_jitter() +
stat_summary(fun.y=mean, geom="point", color="red") +
stat_summary(fun.data=n_fun, geom="text") +
xlab("") +
ylab("N sequences / clone") +
facet_grid(cols=vars(PATIENT), rows = vars(POPULATION), drop=T)
ggsave(plot = pclones2, filename = paste(patdir_number, "/Clones_sequences_population_", clones$PATIENT[1], "_", clones$TREATMENT[1],"grid.svg",sep = ""), device="svg", width = 25, height = 25, units="cm")
ggsave(plot = pclones2, filename = paste(patdir_number, "/Clones_sequences_population_", clones$PATIENT[1], "_", clones$TREATMENT[1],"grid.png",sep = ""), device="png", width = 25, height = 25, units="cm")

clone_N <- count(clones, vars="SAMPLE_POP")
clone_N$vars <- NULL

clone_mean <- tapply(clones$SEQ_COUNT, clones$SAMPLE_POP, mean)
clone_mean <- data.frame(as.list(clone_mean))
clone_mean <- t(clone_mean)
clone_mean <- as.data.frame(clone_mean)
colnames(clone_mean) <- c("mean")
clone_mean$sample <- row.names(clone_mean)
row.names(clone_mean) <- NULL

clone_median <- tapply(clones$SEQ_COUNT, clones$SAMPLE_POP, median)
clone_median <- data.frame(as.list(clone_median))
clone_median <- t(clone_median)
clone_median <- as.data.frame(clone_median)
colnames(clone_median) <- c("median")

clone_stats <- merge(x=clone_mean, y=clone_median, by.x = "sample", by.y = "row.names", all.x=T)
clone_number_stats <- merge(x = clone_N, y=clone_stats, by.x = "SAMPLE_POP", by.y = "sample")
colnames(clone_number_stats) <- c("Sample", "Number of clones", "Mean number of sequences per clone", "Median number of sequences per clone")
write.table(clone_number_stats, file =paste(patdir_number,"/Clones_number_stats_population_", clones$PATIENT[1],"_",clones$TREATMENT[1], ".tsv", sep=""), sep = "\t", row.names = F, quote = F)


##################
# Clonal lineages
#################

# save clonal table
countclones <- countClones(df_pat,clone="CLONE", copy="DUPCOUNT")
write.table(countclones, paste(patdir_lineage, "/", "Clones_table_patient_", df_pat$SOURCE[1],".tsv", sep=""), quote=F, sep="\t", row.names = F)

# Restrict clonal tree size
clones <- subset(countclones, SEQ_COUNT >= 2 & SEQ_COUNT < 1000)

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

    #dnapars_exec <- "/opt/conda/envs/nf-core-bcellmagic-1.1.1dev/bin/dnapars"
    dnapars_exec <- "dnapars"
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

    # Plot
    svg(filename = paste(patdir_lineage_trees,"/Clone_tree_", clone@data$SOURCE[1], "_clone_id_", clone_id, ".svg", sep=""))
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
        vertex.label.color="black", vertex.size=(vsize/20 + 6))
    legend("topleft", c("Germline", "Inferred", "Sample"), 
        fill=c("black", "white", "steelblue"), cex=0.75)
    dev.off()
    
    #Save graph in graphML format
    write_graph(graph, file=paste(patdir_lineage_graphml, "/Graph_", clone@data$SOURCE[1],  "_clone_id_", clone_id, ".txt", sep=""), format = c("graphml"))

}

for (clone_id in clones$CLONE){
    save_graph(df_pat, clone_id)
#    tryCatch(withCallingHandlers(save_graph(df_pat, clone_id), 
#                    error=function(e) {print(paste0("Skipping clone due to problem:", clone_id))},
#                    warning=function(w) {print(paste0("Warning for clone:", clone_id))
#                                invokeRestart("muffleWarning")}), 
#            error = function(e) { print(paste0("Processed clone:", clone_id)) })

}

