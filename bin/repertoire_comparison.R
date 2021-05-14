#!/usr/bin/env Rscript

library(alakazam)
library(shazam)
library(stringr)
library(dplyr)

theme_set(theme_bw(base_family = "ArialMT") + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))



#datadir <- "Data"
datadir <- "."
outdir <- "repertoire_comparison"

# setwd to results folder (containing alakazam, shazam, etc. folders)

### Read all the tables as produced by the pipeline in the current folder and joins them together in the df_all dataframe

all_files <- system(paste0("find '",datadir,"' -name '*germ-pass.tsv'"), intern=T)

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

# Set number of bootrstraps
nboot = 200

# Generate one big dataframe from all patient dataframes
df_all = data.frame()
for (file in all_files){
  fname = file
  print(fname)
  
  df_pat <- read.csv(fname, sep="\t")
  
  df_all <- rbind(df_all, df_pat)
  
}
write.table(df_all, "all_data.tsv", sep = "\t", quote=F, row.names = F, col.names = T)

# Remove underscores in these columns
df_all$treatment <- sapply(df_all$treatment, function(x) str_replace(as.character(x), "_", ""))
df_all$source <- sapply(df_all$source, function(x) str_replace(as.character(x), "_", ""))
df_all$extract_time <- sapply(df_all$extract_time, function(x) str_replace(as.character(x), "_", ""))
df_all$population <- sapply(df_all$population, function(x) str_replace(as.character(x), "_", ""))

# Annotate sample and samplepop (sample + population) by add ing all the conditions
df_all$sample <- as.factor(paste(df_all$treatment, df_all$extract_time, df_all$source, sep="_"))
df_all$sample_pop <- as.factor(paste(df_all$treatment, df_all$extract_time, df_all$source, df_all$population, sep="_"))

##########################
## ABUNDANCE PER PATIENT
#########################

# Per patient
abund <- estimateAbundance(df_all, clone="clone_id", group = "sample", ci=0.95, nboot=nboot)
abund@abundance$treatment <- sapply(abund@abundance$sample, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@abundance$time_point <- sapply(abund@abundance$sample, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@abundance$patient <- sapply(abund@abundance$sample, function(x) unlist(strsplit(as.character(x), "_"))[3])

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")

p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                                  group = sample)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = patient), alpha = 0.4) + 
  geom_line(aes(color = patient)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(cols = vars(treatment), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_patient.pdf"), device="pdf", width = 25, height = 10, units="cm")
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_patient.png"), device="png", width = 25, height = 10, units="cm")

write.table(abund@abundance, file = paste0(abundance_dir, "/Clonal_abundance_data_patient.tsv"), sep="\t", quote = F, row.names = F)

########################
## DIVERSITY PER PATIENT
#######################

print("Diversity calculation")
# Plotting sample diversity per patient
sample_div <- alphaDiversity(abund, group="sample", min_q=0, max_q=4, step_q=0.05,
                             ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n[1], ")")

sample_div@diversity$treatment <- sapply(sample_div@diversity$sample, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@diversity$time_point <- sapply(sample_div@diversity$sample, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@diversity$patient <- sapply(sample_div@diversity$sample, function(x) unlist(strsplit(as.character(x), "_"))[3])

p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                       group = sample)) + 
  geom_ribbon(aes(ymin = d_lower, 
                  ymax = d_upper, 
                  fill = patient), alpha = 0.4) +
  geom_line(aes(color = patient)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_grid(cols=vars(treatment))
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.png"), device="png", width = 25, height = 10, units="cm")
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.pdf"), device="pdf", width = 25, height = 10, units="cm")

# DIVERSITY AT Q=1

sample_div_q1 <- sample_div@diversity[which(sample_div@diversity$q == 1),]

sample_main <- paste0("Sample diversity at Q=1 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=patient)) + 
  geom_point(position=dodge, stat="identity") +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(treatment), drop=T, space="free", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.png"), device="png", 
       width = 20, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.pdf"), device="pdf", 
       width = 20, height = 10, units="cm")

# DIVERSITY AT Q=0

sample_div_q0 <- sample_div@diversity[which(sample_div@diversity$q == 0),]

sample_main <- paste0("Sample diversity at Q=0 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=patient)) + 
  geom_point(position=dodge, stat="identity") +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(treatment), drop=T, space="free", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient.png"), device="png", 
       width = 20, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient.pdf"), device="pdf", 
       width = 20, height = 10, units="cm")

sample_test <- sample_div@tests
write.table(sample_test, file = paste0(diversity_dir, "/Diversity_tests_data_patient.tsv"), 
            sep="\t", quote = F, row.names = F)

##########################
# ABUNDANCE PER POPULATION
##########################


abund <- estimateAbundance(df_all, clone="clone_id", group = "sample_pop", ci=0.95, nboot=nboot)
abund@abundance$treatment <- sapply(abund@abundance$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@abundance$time_point <- sapply(abund@abundance$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@abundance$patient <- sapply(abund@abundance$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[3])
abund@abundance$population <- sapply(abund@abundance$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[4])

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")
p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                                  group = sample_pop)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = population), alpha = 0.4) + 
  geom_line(aes(color = population)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(vars(patient), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_patient_population.pdf"), device="pdf", 
       width = 30, height = 20, units="cm")
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_patient_population.png"), device="png", 
       width = 30, height = 20, units="cm")


write.table(abund@abundance, file = paste0(abundance_dir, "/Clonal_abundance_data_patient_population.tsv"), sep="\t", 
            quote = F, row.names = F)

##########################
# ABUNDANCE PER STATE
##########################

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")
p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                                  group = sample_pop)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = patient), alpha = 0.4) + 
  geom_line(aes(color = patient)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(rows=vars(population), cols=vars(treatment), drop=T, scales="free") +
  #facet_wrap(vars(population,treatment), drop = T, nrow=4, ncol=4) + 
  theme(strip.text.y = element_text(size = 5))
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_population_state.svg"), device="svg", 
       width = 15, height = 20, units="cm")
ggsave(plot=p1, filename = paste0(abundance_dir,"/Clonal_abundance_population_state.png"), device="png",
       width = 15, height = 20, units="cm")

###########################
# DIVERSITY PER POPULATION
###########################

# Plotting sample diversity for all populations
print("Diversity calculation population")

sample_div <- alphaDiversity(abund, group="sample_pop", min_q=0, max_q=4, step_q=0.05,
                             ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n[1], ")")

sample_div@diversity$treatment <- sapply(sample_div@diversity$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@diversity$time_point <- sapply(sample_div@diversity$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@diversity$patient <- sapply(sample_div@diversity$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[3])
sample_div@diversity$population <- sapply(sample_div@diversity$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[4])

p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                       group = sample_pop)) + 
  geom_ribbon(aes(ymin = d_lower, 
                  ymax = d_upper, fill = population), alpha = 0.4) +
  geom_line(aes(color = population)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_wrap(vars(patient), scales = "free_x")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.svg"), device="svg", 
       width = 25, height = 20, units="cm")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.pdf"), device="pdf", 
       width = 25, height = 20, units="cm")

###########################
# DIVERSITY PER STATE
###########################
p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                       group = sample_pop)) + 
  geom_ribbon(aes(ymin = d_lower, 
                  ymax = d_upper, fill = patient), alpha = 0.4) +
  geom_line(aes(color = patient)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) +
  theme(strip.text.y = element_text(size = 5)) +   
  facet_grid(rows=vars(population), cols=vars(treatment), drop=T, scales="free_x")
#facet_wrap(vars(population,treatment), scales = "free_x")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_population_state.svg"), device="svg", 
       width = 30, height = 20, units="cm")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_population_state.pdf"), device="pdf", 
       width = 30, height = 20, units="cm")

# Tests sample diversity for significance
print("Diversity calculation tests population")

# DIVERSITY AT Q=1

sample_div_q1 <- sample_div@diversity[which(sample_div@diversity$q == 1),]

paste0("Sample diversity at Q=1 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=population, fill=population)) + 
  geom_point(position=dodge, stat="identity", aes(fill=population)) +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(patient)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.svg"), device="svg", 
       width = 30, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.pdf"), device="pdf", 
       width = 30, height = 10, units="cm")

#Diversity Q1 by state
dodge <- position_dodge(width = 0.9)
g2 <- ggplot(sample_div_q1, aes(y=d, x=treatment)) + 
  geom_jitter(width=0.05,stat="identity", aes(color=patient)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, size = 0.3) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  theme(axis.text.x = element_text(angle = 49, hjust = 1), legend.position = "right")  +
  facet_wrap(vars(population), scales="free", drop = T)

ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q1_population_state.svg"), device="svg", 
       width = 30, height = 20, units="cm")
ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q1_population_state.pdf"), device="pdf", 
       width = 30, height = 20, units="cm")

# DIVERSITY AT Q=0

sample_div_q0 <- sample_div@diversity[which(sample_div@diversity$q == 0),]

paste0("Sample diversity at Q=0 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=population, fill=population)) + 
  geom_point(position=dodge, stat="identity", aes(fill=population)) +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(patient)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient_population.svg"), device="svg", 
       width = 30, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient_population.pdf"), device="pdf", 
       width = 30, height = 10, units="cm")

dodge <- position_dodge(width = 0.9)
g2 <- ggplot(sample_div_q0, aes(y=d, x=treatment)) + 
  geom_jitter(width=0.05, stat="identity", aes(color=patient)) +
  # geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_wrap(vars(population), scales="free_x", drop = T) +
  # facet_grid(cols=vars(treatment), rows=vars(population)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, size = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") 

ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q0_population_state.svg"), device="svg", 
       width = 30, height = 15, units="cm")
ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q0_population_state.pdf"), device="pdf", 
       width = 30, height = 15, units="cm")


sample_test <- sample_div@tests
write.table(sample_test, file = paste0(diversity_dir, "/Diversity_tests_data_patient_population.tsv"), 
            sep="\t", quote = F, row.names = F)


############
## ISOTYPES
############

print("Isotypes calculation")
# Plotting Isotype percentages per patient
df_all$isotype <- df_all$c_primer

res <- df_all %>% group_by(isotype,sample,source,treatment,extract_time) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(source),])
res_sample <- df_all %>% group_by(sample) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "sample", by.y = "sample")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=isotype, y=Freq, x=isotype)) +
  geom_bar(position = "dodge", stat="identity") +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  facet_wrap(vars(source), scales="free_x", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.svg"), device = "svg", 
       width = 18, height = 15, units = "cm")
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.pdf"), device = "pdf", 
       width = 18, height = 15, units = "cm")

write.table(freqs, file = paste0(isotype_dir,"/Isotype_frequencies_data.tsv"), sep="\t", quote=F, row.names = F)

# Plotting Isotype percentages per population
print("Isotypes calculation per population")
res <- df_all %>% group_by(isotype,sample_pop,source,treatment,extract_time, population) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(source),])
res_sample <- df_all %>% group_by(sample_pop) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "sample_pop", by.y = "sample_pop")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=isotype, y=Freq, x=isotype)) +
  geom_bar(position = "dodge", stat="identity") +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  facet_grid(cols=vars(source), rows=vars(population)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.svg"), device = "svg", 
       width = 25, height = 20, units = "cm")
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.pdf"), device = "pdf", 
       width = 25, height = 20, units = "cm")

write.table(freqs, file = paste0(isotype_dir, "/Isotype_frequencies_population_data.tsv"), sep="\t", quote = F, row.names = F)

#Plotting Isotype percentages per state
g4 <- ggplot(freqs, aes(y=Freq, x=treatment)) +
  geom_point(position = "dodge", stat="identity", aes(color=source)) +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, size = 0.3) +
  facet_grid(cols =vars(isotype), rows=vars(population),scales="free", drop = T)+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population_state.svg"), device = "svg", 
       width = 25, height = 12, units = "cm")
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population_state.pdf"), device = "pdf", 
       width = 25, height = 12, units = "cm")

##################
# V FAMILY USAGE              
##################

# Quantify V family usage by patient
print("V family usage calculation")
family <- countGenes(df_all, gene="v_call", groups="sample", 
                     mode="family")
family$treatment <- sapply(family$sample, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$time_point <- sapply(family$sample, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$patient <- sapply(family$sample, function(x) unlist(strsplit(as.character(x), "_"))[3])

g2 <- ggplot(family, aes(x=gene, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=gene), size=4, alpha=0.8) +
  ggtitle("V Gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_wrap(vars(patient), scales="free_x")
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.svg"), plot = g2, width = 18, height = 15, units = "cm")
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.pdf"), plot = g2, width = 18, height = 15, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data.tsv"), sep = "\t", quote = F, row.names = F)

# Quantify V family usage by population
print("V family usage calculation population")
family <- countGenes(df_all, gene="v_call", groups="sample_pop", 
                     mode="family")
family$treatment <- sapply(family$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$time_point <- sapply(family$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$patient <- sapply(family$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[3])
family$population <- sapply(family$sample_pop, function(x) unlist(strsplit(as.character(x), "_"))[4])

g2 <- ggplot(family, aes(x=gene, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=gene), size=3, alpha=0.8) +
  ggtitle("V gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(patient), rows=vars(population))
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.svg"), plot = g2, 
       width = 30, height = 20, units = "cm")
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.pdf"), plot = g2, 
       width = 30, height = 20, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data_population.tsv"), sep = "\t", 
            quote = F, row.names = F)

# Quantify V family usage by state
g2 <- ggplot(family, aes(x=treatment, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=patient), size=3, alpha=0.8) +
  ggtitle("V gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(population), rows=vars(gene),scales="free", drop = T) 
#facet_grid(cols=vars(patient,treatment), rows=vars(population))
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population_state.svg"), plot = g2, 
       width = 30, height = 20, units = "cm")
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population_state.pdf"), plot = g2, 
       width = 30, height = 20, units = "cm")

