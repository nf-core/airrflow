#!/usr/bin/env Rscript

library("alakazam")
library("ggplot2")
library("data.table")
library("dplyr")
library("tigger")
library("shazam")
library("igraph")
library("svglite")
library("extrafont")

#extrafont::font_import(prompt = FALSE)
#extrafont::loadfonts()

theme_set(theme_bw(base_family = "ArialMT") + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))



datadir <- "."
outdir <- "repertoire_analysis"

# setwd to results folder (containing alakazam, shazam, etc. folders)

### Read all the tables as produced by the pipeline in the current folder and joins them together in the df_all dataframe

all_files <- system(paste0("find '",datadir,"' -name '*.tab'"), intern=T)

dir.create(outdir)
diversity_dir <- paste(outdir, "Diversity", sep="/")
abundance_dir <- paste(outdir, "Abundance", sep="/")
isotype_dir <- paste(outdir, "Isotype", sep="/")
vfamily_dir <- paste(outdir, "V_family", sep="/")
mutation_dir <- paste(outdir, "Mutational_load", sep="/")
dir.create(diversity_dir)
dir.create(abundance_dir)
dir.create(isotype_dir)
dir.create(vfamily_dir)
dir.create(mutation_dir)

nboot = 1000

df_all = data.frame()
for (file in all_files){
   fname = file
   print(fname)

   df_pat <- read.csv(fname, sep="\t")
   
   df_all <- rbind(df_all, df_pat)

}


# Annotate sample and samplepop (sample + population) by adding all the conditions
df_all$SAMPLE <- as.factor(paste(df_all$TREATMENT, df_all$EXTRACT_TIME, df_all$SOURCE, sep="_"))
df_all$SAMPLE_POP <- as.factor(paste(df_all$TREATMENT, df_all$EXTRACT_TIME, df_all$SOURCE, df_all$POPULATION, sep="_"))

###############
## DIVERSITY
###############

print("Diversity calculation")
# Plotting sample diversity per patient
sample_div <- rarefyDiversity(df_all, "SAMPLE", min_q=0, max_q=4, step_q=0.05,
                              ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")

sample_div@data$TREATMENT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@data$TIME_POINT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@data$PATIENT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])

p2 <- ggplot(sample_div@data, aes(x = Q, y = D, 
                                  group = GROUP)) + 
  geom_ribbon(aes(ymin = D_LOWER, 
                  ymax = D_UPPER, 
                  fill = TIME_POINT), alpha = 0.4) +
  geom_line(aes(color = TIME_POINT, linetype=PATIENT)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_grid(cols=vars(TREATMENT))
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.svg"), device="svg", width = 25, height = 7, units="cm")
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.pdf"), device="pdf", width = 25, height = 7, units="cm")


# Tests sample diversity for significance per patient
# print("Diversity calculation tests")

# sample_test <- testDiversity(df_all, 1, "SAMPLE", nboot=nboot)

# sample_test@summary$TREATMENT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
# sample_test@summary$TIME_POINT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
# sample_test@summary$PATIENT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])

# dodge <- position_dodge(width = 0.9)
# g1 <- ggplot(sample_test@summary, aes(y=MEAN, x=PATIENT, fill=TIME_POINT)) + 
#   geom_bar(position=dodge, stat="identity") +
#   geom_errorbar(aes(ymin=MEAN-SD, ymax=MEAN+SD), width = .2, position=dodge) +
#   xlab("") + ylab("Diversity (q=1)") +
#   ggtitle(sample_main) +
#   facet_grid(cols=vars(TREATMENT), drop=T, space="free", scales = "free") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.svg"), device="svg", 
# width = 25, height = 10, units="cm")
# ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.pdf"), device="pdf", 
# width = 25, height = 10, units="cm")

# write.table(sample_test@summary, file = paste0(diversity_dir, "/Diversity_Q1_data_patient.tsv"), 
# sep="\t", quote = F, row.names = F)


# Plotting sample diversity for all populations
print("Diversity calculation population")

sample_div <- rarefyDiversity(df_all, "SAMPLE_POP", min_q=0, max_q=4, step_q=0.05,
                              ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")

sample_div@data$TREATMENT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@data$TIME_POINT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@data$PATIENT <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])
sample_div@data$POPULATION <- sapply(sample_div@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[4])

p2 <- ggplot(sample_div@data, aes(x = Q, y = D, 
                                   group = GROUP)) + 
  geom_ribbon(aes(ymin = D_LOWER, 
            ymax = D_UPPER, fill = TIME_POINT), alpha = 0.4) +
  geom_line(aes(color = TIME_POINT, linetype = PATIENT)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_grid(cols=vars(TREATMENT), rows=vars(POPULATION))
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.svg"), device="svg", 
width = 25, height = 20, units="cm")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.pdf"), device="pdf", 
width = 25, height = 20, units="cm")

# Tests sample diversity for significance
# print("Diversity calculation tests population")
# sample_test <- testDiversity(df_all, 1, "SAMPLE_POP", nboot=nboot)

# sample_test@summary$TREATMENT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
# sample_test@summary$TIME_POINT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
# sample_test@summary$PATIENT <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])
# sample_test@summary$POPULATION <- sapply(sample_test@summary$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[4])

# dodge <- position_dodge(width = 0.9)
# g1 <- ggplot(sample_test@summary, aes(y=MEAN, x=PATIENT, fill=TIME_POINT)) + 
#   geom_bar(position=dodge, stat="identity") +
#   geom_errorbar(aes(ymin=MEAN-SD, ymax=MEAN+SD), width = .2, position=dodge) +
#   xlab("") + ylab("Diversity (q=1)") +
#   ggtitle(sample_main) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_grid(cols=vars(TREATMENT), rows=vars(POPULATION), drop = T, scales = "free")
# ggsave(plot=g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.svg"), device="svg", 
# width = 25, height = 10, units="cm")
# ggsave(plot=g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.pdf"), device="pdf", 
# width = 25, height = 10, units="cm")

# write.table(sample_test@summary, file = paste0(diversity_dir, "/Diversity_Q1_data_patient_population.tsv"), 
#   sep="\t", quote = F, row.names = F)

####################
## ABUNDANCE
####################

# Per patient
abund <- estimateAbundance(df_all, group = "SAMPLE", ci=0.95, nboot=nboot)
abund@data$TREATMENT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@data$TIME_POINT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@data$PATIENT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])

abund_main <- paste0("Rank abundance")

p1 <- ggplot(abund@data, aes(x = RANK, y = P, 
                                   group = GROUP)) + 
  geom_ribbon(aes(ymin = LOWER, 
                         ymax = UPPER, fill = TIME_POINT), alpha = 0.4) + 
  geom_line(aes(color = TIME_POINT, linetype = PATIENT)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(cols = vars(TREATMENT), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_patient.svg"), device="svg", width = 25, height = 10, units="cm")

write.table(abund@data, file = paste0(abundance_dir, "/Abundance_data_patient.tsv"), sep="\t", quote = F, row.names = F)

# Per population
abund <- estimateAbundance(df_all, group = "SAMPLE_POP", ci=0.95, nboot=nboot)
abund@data$TREATMENT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@data$TIME_POINT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@data$PATIENT <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[3])
abund@data$POPULATION <- sapply(abund@data$GROUP, function(x) unlist(strsplit(as.character(x), "_"))[4])

abund_main <- paste0("Rank abundance")

p1 <- ggplot(abund@data, aes(x = RANK, y = P, 
                             group = GROUP)) + 
  geom_ribbon(aes(ymin = LOWER, 
                  ymax = UPPER, fill = TIME_POINT), alpha = 0.4) + 
  geom_line(aes(color = TIME_POINT, linetype = PATIENT)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(cols = vars(TREATMENT), rows = vars(POPULATION), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_patient_population.svg"), device="svg", 
width = 25, height = 25, units="cm")

write.table(abund@data, file = paste0(abundance_dir, "/Abundance_data_patient_population.tsv"), sep="\t", 
quote = F, row.names = F)

############
## ISOTYPES
############

print("Isotypes calculation")
# Plotting Isotype percentages per patient
df_all$ISOTYPE <- df_all$C_PRIMER

res <- df_all %>% group_by(ISOTYPE,SAMPLE,SOURCE,TREATMENT,EXTRACT_TIME) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(SOURCE),])
res_sample <- df_all %>% group_by(SAMPLE) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "SAMPLE", by.y = "SAMPLE")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=EXTRACT_TIME, y=Freq, x=ISOTYPE)) +
  geom_bar(position = "dodge", stat="identity") +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  facet_grid(cols=vars(TREATMENT, SOURCE), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.svg"), device = "svg", 
  width = 25, height = 7, units = "cm")
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.pdf"), device = "pdf", 
  width = 25, height = 7, units = "cm")

write.table(freqs, file = paste0(isotype_dir,"/Isotype_frequencies_data.tsv"), sep="\t", quote=F, row.names = F)

# Plotting Isotype percentages per population
print("Isotypes calculation per population")
res <- df_all %>% group_by(ISOTYPE,SAMPLE_POP,SOURCE,TREATMENT,EXTRACT_TIME, POPULATION) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(SOURCE),])
res_sample <- df_all %>% group_by(SAMPLE_POP) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "SAMPLE_POP", by.y = "SAMPLE_POP")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=EXTRACT_TIME, y=Freq, x=ISOTYPE)) +
 geom_bar(position = "dodge", stat="identity") +
 xlab("") + 
 ylab("Frequency") +
 ggtitle("Isotype frequency") +
 facet_grid(cols=vars(TREATMENT, SOURCE), rows=vars(POPULATION)) +
 theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.svg"), device = "svg", 
  width = 25, height = 20, units = "cm")
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.pdf"), device = "pdf", 
  width = 25, height = 20, units = "cm")

write.table(freqs, file = paste0(isotype_dir, "/Isotype_frequencies_population_data.tsv"), sep="\t", quote = F, row.names = F)

##################
# V FAMILY USAGE              
##################

# Quantify V family usage by patient
print("V family usage calculation")
family <- countGenes(df_all, gene="V_CALL", groups="SAMPLE", 
                     mode="family")
family$TREATMENT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$TIME_POINT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$PATIENT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[3])

g2 <- ggplot(family, aes(x=GENE, y=SEQ_FREQ)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=TIME_POINT), size=5, alpha=0.8) +
  ggtitle("V Gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(TREATMENT, PATIENT))
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.svg"), plot = g2, width = 25, height = 10, units = "cm")
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.pdf"), plot = g2, width = 25, height = 10, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data.tsv"), sep = "\t", quote = F, row.names = F)

# Quantify V family usage by population
print("V family usage calculation population")
family <- countGenes(df_all, gene="V_CALL", groups="SAMPLE_POP", 
                     mode="family")
family$TREATMENT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$TIME_POINT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$PATIENT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[3])
family$POPULATION <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[4])

g2 <- ggplot(family, aes(x=GENE, y=SEQ_FREQ)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=TIME_POINT), size=5, alpha=0.8) +
  ggtitle("V gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(PATIENT, TREATMENT), rows=vars(POPULATION))
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.svg"), plot = g2, 
  width = 25, height = 25, units = "cm")
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.pdf"), plot = g2, 
  width = 25, height = 25, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data_population.tsv"), sep = "\t", 
  quote = F, row.names = F)

#######################
# MUTATIONAL FREQUENCY
#######################

# Quantify mutational load by patient
print("Mutational load calculation")
# Calculate mutation counts
df_all_mut_counts <- observedMutations(df_all, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     regionDefinition=NULL,
                                     frequency=FALSE,
                                     nproc=4) 
# Calculate mutation frequencies
df_all_mut_freq <- observedMutations(df_all_mut_counts, sequenceColumn = "SEQUENCE_IMGT",
                                     germlineColumn = "GERMLINE_IMGT_D_MASK",
                                     regionDefinition = NULL,
                                     frequency = TRUE,
                                     nproc = 4)

mut_counts_freqs <- select(df_all_mut_freq, SEQUENCE_ID, starts_with("MU_FREQ_"), starts_with("MU_COUNT_"))

res_mut <- mut_counts_freqs %>% group_by(SAMPLE,SOURCE,TREATMENT,EXTRACT_TIME) %>% 
  dplyr::summarise(MUTATION_MEAN_COUNT=mean(MU_COUNT), 
  MUTATION_MEDIAN_COUNT=median(MU_COUNT), MUTATION_SD_COUNT=sd(MU_COUNT), MUTATION_MEAN_FREQ=mean(MU_FREQ), MUTATION_MEDIAN_FREQ=median(MU_FREQ), MUTATION_SD_FREQ=mean(MU_FREQ), N_SEQS = n())
write.table(res_mut, file = paste0(mutation_dir,"/Mutation_stats_patient.tsv"), sep="\t", col.names = T, row.names = F, quote = F)

res_mut_pop <- mut_counts_freqs %>% group_by(SAMPLE_POP,SOURCE,TREATMENT,EXTRACT_TIME, POPULATION) %>% 
  dplyr::summarise(MUTATION_MEAN_COUNT=mean(MU_COUNT), MUTATION_MEDIAN_COUNT=median(MU_COUNT), MUTATION_SD_COUNT=sd(MU_COUNT), 
  MUTATION_MEAN_FREQ=mean(MU_FREQ), MUTATION_MEDIAN_FREQ=median(MU_FREQ), MUTATION_SD_FREQ=mean(MU_FREQ), N_SEQS = n())
write.table(res_mut_pop, file = paste0(mutation_dir,"/Mutation_stats_patient_population.tsv"), sep="\t", col.names = T, row.names = F, quote = F)


plot_mut_num <- ggplot(mut_counts_freqs, aes(fill=EXTRACT_TIME, y=MU_COUNT, x=SAMPLE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Counts") +
  ggtitle("Mutation Counts") +
  facet_grid(cols=vars(TREATMENT, SOURCE), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient.svg"), device = "svg", width = 25, height = 7, units = "cm")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient.png"), device = "png", width = 25, height = 7, units = "cm")


plot_mut_freq <- ggplot(mut_counts_freqs, aes(fill=EXTRACT_TIME, y=MU_FREQ, x=SAMPLE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Frequency") +
  ggtitle("Mutation Frequency") +
  facet_grid(cols=vars(TREATMENT, SOURCE), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient.svg"), device = "svg", width = 25, height = 7, units = "cm")
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient.png"), device = "png", width = 25, height = 7, units = "cm")


plot_mut_num <- ggplot(mut_counts_freqs, aes(fill=EXTRACT_TIME, y=MU_COUNT, x=SAMPLE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Counts") +
  ggtitle("Mutation Counts per Population") +
  facet_grid(cols=vars(TREATMENT, SOURCE), rows=vars(POPULATION), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient_population.svg"), device = "svg", width = 25, height = 20, units = "cm")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient_population.png"), device = "png", width = 25, height = 20, units = "cm")


plot_mut_freq <- ggplot(mut_counts_freqs, aes(fill=EXTRACT_TIME, y=MU_FREQ, x=SAMPLE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Frequency") +
  ggtitle("Mutation Frequency per Population") +
  facet_grid(cols=vars(TREATMENT, SOURCE), rows=vars(POPULATION), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient_population.svg"), device = "svg", width = 25, height = 20, units = "cm")
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient_population.png"), device = "png", width = 25, height = 20, units = "cm")
