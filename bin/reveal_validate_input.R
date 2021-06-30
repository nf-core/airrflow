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
# Valid species. Use lower case. 
valid_species <- list("human"=c("human", "homo sapiens", "h. sapiens", "hs"),
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
miairr_metadata <- read.csv(opt$MIAIRR, sep="\t", stringsAsFactors = F)
mandatory <- miairr_metadata[['Mandatory.BioSample.attribute']] == TRUE
mandatory_fields <- setdiff(c(miairr_metadata[['AIRR.Formats.WG.field.name']][mandatory], "pcr_target_locus", "single_cell", "species"), "organism")

if (!all(mandatory_fields %in% colnames(input))) {
   missing_fields <- mandatory_fields[mandatory_fields %in% colnames(input) == FALSE]
   if ((missing_fields == "species") & ("organism" %in% colnames(input))) {
         message("Missing 'species' field, using 'organism' field instead.")
         input[['species']] <- input[['organism']]
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

# `species` must exist, then translate to name format 
# suitable for igblast and downstream analysis
get_valid_species <- function(dictionary, words) {
   dictionary <-  bind_rows(lapply(dictionary, function(x) {data.frame("species"=x)} ), .id="valid_species")
   sapply(words, function(word){
      valid_species <- dictionary %>%
         filter(species == tolower(word)) %>%
         distinct(valid_species)
      if (nrow(valid_species)>0){
         if(length(col)>1) {
            warning("Multiple `species` mappings found for ",word,": ", paste(as.character(valid_species), collapse=","))
            NA
         } else {
            as.character(valid_species)
         }
      } else {
         warning("No `species` mapping available for: `",word,"`.")
         NA
      } 
   })
}

species <- get_valid_species(valid_species, input[['species']])
input[,'species'] <- species
not_valid_species <- species[is.na(species)]

input[,'valid_species'] <- TRUE
if (length(not_valid_species)>0) {
   input[is.na(species),'valid_species'] <- FALSE  
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

input$id <- paste0("input_id_", 1:nrow(input))

# Validate single_Cell field. Should be true or false
input$valid_single_cell <- input[['single_cell']] %in% c(T,F,"true", "false", "TRUE","FALSE")
is_single_cell <- input[['single_cell']] %in% c(T,"true", "TRUE") 
if (sum(is_single_cell)>0) {
   input[['single_cell']][is_single_cell] <- "true"
}
if (sum(!is_single_cell)>0) {
   input[['single_cell']][!is_single_cell] <- "false"
}

# Validate pcr_target_locus
input[['valid_pcr_target_locus']]  <- tolower(input[['pcr_target_locus']]) %in% tolower(c("IGH", "IGI", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG", "IG", "TR"))
input[['locus']] <- tolower(substr(input[['pcr_target_locus']],1,2))

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

check_fields <- c("valid_filename", "valid_species", "valid_collapseby", "valid_cloneby", "valid_single_cell", "valid_pcr_target_locus")
not_valid_rows <- which(rowSums(input[, check_fields, drop=FALSE]) < length(check_fields))

if (length(not_valid_rows) > 0) {
   write.table(input[not_valid_rows,], file=paste0(opt$OUTPUT,"_not-valid.tsv"), sep="\t",
               quote=FALSE, row.names = FALSE)
   stop("Please review ", opt$INPUT,". Invalid information found. Inspect ", paste0(opt$OUTPUT,".html"), ".")
}
