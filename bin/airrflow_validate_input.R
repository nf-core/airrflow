#!/usr/bin/env Rscript
# 
# Validate --input. Check: 
# - filenames are .fasta, .fa or .tsv
# - column names are MiAIRR compliant
#
# Arguments:
#   --input   Tabulated data
#   --collapseby Names of the metadata columns to be used to collapse
#              duplicated sequences. Default: sample_id
#   --cloneby Names of the metadata columns that need to be used to group
#             files to identify clonal groups. Default: subject_id
#   --output  Output name. Default:validate_input
#   -h  Display help.
# Example: ./validate_input.R --input ../../test-datasets/metadata.tsv --collapseby sample_id --cloneby subject_id

# TODO: should validate metadata using AIRR schema, and properly tested in the
# framework, so that we incorporate any updates in the standard
# https://github.com/airr-community/airr-standards/blob/v1.3.1/specs/airr-schema.yaml
# https://github.com/airr-community/airr-standards/blob/v1.3.1/NCBI_implementation/mapping_MiAIRR_BioSample.tsv
# metadata data.frame
# miairr_mapping. Local version of https://github.com/airr-community/airr-standards/blob/v1.3.1/NCBI_implementation/mapping_MiAIRR_BioSample.tsv

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))

# Define commmandline arguments
opt_list <- list(make_option(c("--input"), dest="INPUT", default=NULL,
                             help="Input file."),
                 make_option(c("--collapseby"), dest="COLLAPSEBY", default='sample_id',
                             help="Grouping fields to collapse duplicated sequences."),                    
                 make_option(c("--cloneby"), dest="CLONEBY", default='subject_id',
                             help="Grouping fields to identify clonally related sequences."),                 
                 make_option(c("--output"), dest="OUTPUT", default="validate_input",
                             help="Output name."),
                 make_option(c("--miairr"), dest="MIAIRR", default=NULL,
                             help="Output name."))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))
opt
getwd()

# Settings
# Valid organisms. Use lower case. 
valid_organisms <- list("human"=c("human", "homo sapiens", "h. sapiens", "hs"),
                        "mouse"=c("mouse", "mus musculus", "mm"))
# Check miairr file
if (!("MIAIRR" %in% names(opt))) {
   stop("You must provide a miairr standard file with --miairr.")
}

# Check input file
if (!("INPUT" %in% names(opt))) {
    stop("You must provide an input file with --input.")
}

# Read metadata file
input <- read.csv(opt$INPUT,sep = "\t", header=TRUE, stringsAsFactors = F)

# Validate metadata file
# TODO: organism not used anymore, standard changed to species
miairr_metadata <- read.csv(opt$MIAIRR, sep="\t", stringsAsFactors = F)
mandatory <- miairr_metadata[['Mandatory.BioSample.attribute']] == TRUE
mandatory_fields <- miairr_metadata[['AIRR.Formats.WG.field.name']][mandatory]

if (!all(mandatory_fields %in% colnames(input))) {
   missing_fields <- mandatory_fields[mandatory_fields %in% colnames(input) == FALSE]
   if ((missing_fields == "organism") & ("species" %in% colnames(input))) {
         message("Missing 'organism' field, using 'species' field instead.")
         input[['organism']] <- input[['species']]
   } else {
      stop("Missing MiAIRR fields: ", paste(missing_fields,collapse=", "))
   }
}

# Add filetype
input$filetype <- sub(".+\\.([^\\.]+)$","\\1", input$filename)

# Validate `filename`
# `filename` must be unique and must contain .fasta, fa or .tsv files
input$valid_filename <- sapply(basename(input[,'filename']), function(fn) {
   fn <- strsplit(fn, "\\.")[[1]]
   fn[length(fn)] %in% c("fasta", "fa", "tsv")
})

dup_filenames <- input[,'filename'][duplicated(input[,'filename'])]
if (length(dup_filenames)>0) {
   input$valid_filename[input[,'filename'] %in% dup_filenames] <- FALSE
}

#TODO validate MiAIRR fields

# `organism` must exist, then translate to name format 
# suitable for igblast and downstream analysis
get_valid_organism <- function(dictionary, words) {
   dictionary <-  bind_rows(lapply(dictionary, function(x) {data.frame("organism"=x)} ), .id="valid_organism")
   sapply(words, function(word){
      valid_organism <- dictionary %>%
         filter(organism == tolower(word)) %>%
         distinct(valid_organism)
      if (nrow(valid_organism)>0){
         if(length(col)>1) {
            warning("Multiple `organism` mappings found for ",word,": ", paste(as.character(valid_organism), collapse=","))
            NA
         } else {
            as.character(valid_organism)
         }
      } else {
         warning("No `organism` mapping available for: `",word,"`.")
         NA
      } 
   })
}

organisms <- get_valid_organism(valid_organisms, input[['organism']])
input[,'organism'] <- organisms
not_valid_organisms <- organisms[is.na(organisms)]

input[,'valid_organism'] <- TRUE
if (length(not_valid_organisms)>0) {
   input[is.na(organisms),'valid_organism'] <- FALSE  
}


## TODO: make a function to do column name validations
# Validate collapseby
collapseby <- strsplit(opt$COLLAPSEBY,",")[[1]]

# cell_id is a rearrangement property, don't check here
# will be validated when rearrangements are loaded
collapseby <- setdiff(collapseby, "cell_id")

unknown_collapseby <- collapseby[collapseby %in% colnames(input) == FALSE]
if (length(unknown_collapseby)>0) {
   stop("Unknown `collapseby` found: ", paste(unknown_collapseby, collapse = ","))
} else {
   input$valid_collapseby <- T
}

input$collapseby_group <- apply(input,1, function(x) {
   paste(x[collapseby],collapse = ",")
})

input <- input %>%
   group_by(collapseby_group) %>%
   mutate(collapseby_size=n()) %>%
   ungroup()

# Validate cloneby
cloneby <- strsplit(opt$CLONEBY,",")[[1]]
unknown_cloneby <- cloneby[cloneby %in% colnames(input) == FALSE]
if (length(unknown_cloneby)>0) {
   stop("Unknown `cloneby` found: ", paste(unknown_cloneby, collapse = ","))
} else {
   input$valid_cloneby <- T
}

input$cloneby_group <- apply(input,1, function(x) {
   paste(x[cloneby],collapse = ",")
})

input <- input %>%
   group_by(cloneby_group) %>%
   mutate(cloneby_size=n()) %>%
   ungroup()

input$input_id <- paste0("input_id_", 1:nrow(input))

# Write validation results to a file
write.table(input, file=paste0(opt$OUTPUT,".tsv"), sep="\t",
            quote=FALSE, row.names = FALSE)

dt <- DT::datatable(input, 
              rownames = FALSE
              ) %>%
   formatStyle(
      'valid_filename',
      target = 'row',
      backgroundColor = styleEqual(c(0, 1), c('lightsalmon', 'white'))
   )

DT::saveWidget(dt, file=paste0(opt$OUTPUT,".html"))

check_fields <- c("valid_filename", "valid_organism", "valid_collapseby", "valid_cloneby")
not_valid_rows <- which(rowSums(input[, check_fields, drop=FALSE]) < length(check_fields))

if (length(not_valid_rows) > 0) {
   write.table(input[not_valid_rows,], file=paste0(opt$OUTPUT,"_not-valid.tsv"), sep="\t",
               quote=FALSE, row.names = FALSE)
   stop("Please review ", opt$INPUT,". Invalid information found. Inspect ", paste0(opt$OUTPUT,".html"), ".")
}
