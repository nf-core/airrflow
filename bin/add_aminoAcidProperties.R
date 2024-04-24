#!/usr/bin/env Rscript
# Written by Justine Pollet


# Add amino acid properties to repertoire:
# Arguments:
#   --repertoire    Tabulated data in AIRR (TSV) format with clonal assignments, germline assignments and metadata.
#   --outname       Filename for the output repertoire
#   -h  Display help.
# Example: ./add_alakazam_AminoAcidProperties.R --repertoire meta-pass.tsv --outname my-repertoire

# Libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--outname"), default=NULL)
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Read repertoire
repertoire <- read.csv(opt$repertoire, sep="\t", header=TRUE, stringsAsFactors = F)

parsed_sequences=c("cdr1","cdr2","cdr3","fwr1","fwr2","fwr3")

for (seq in 1:length(parsed_sequences)) {

    # Compute amino acid properties for current parsed sequences
    aAp <- aminoAcidProperties(repertoire, seq = parsed_sequences[seq], nt = TRUE, trim = FALSE) %>% 
               select(starts_with(parsed_sequences[seq])) %>%
               select(-parsed_sequences[seq])

    if (seq==1){

        # Add amino acid properties for current parsed sequences into repertoire
        idx_parsed_sequences <- which(colnames(repertoire) == parsed_sequences[seq])
        anno_repertoire <- cbind(repertoire[, 1:idx_parsed_sequences], aAp, repertoire[, (idx_parsed_sequences + 1):ncol(repertoire)])
    }
    if(seq>1){
        # Add amino acid properties for current parsed sequences into repertoire
        idx_parsed_sequences <- which(colnames(anno_repertoire) == parsed_sequences[seq])
        anno_repertoire <- cbind(anno_repertoire[, 1:idx_parsed_sequences], aAp, anno_repertoire[, (idx_parsed_sequences + 1):ncol(anno_repertoire)])
    }
}

# save repertoire table with amino acid properties fields
write.table(anno_repertoire, opt$outname, quote=F, sep="\t", row.names = F)
