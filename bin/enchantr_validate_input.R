#' Create an Immcantation input validation project
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Input Validation
#' to create the skeleton of an Immcantation Input Validation project
#' @param  path path to the directory where the project will be created
validate_input_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "input_validation_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
}

#' @export
validate_input <- function(input, miairr, collapseby, cloneby, reassign=TRUE) {


    # Read metadata file
    input <- read.csv(input,sep = "\t", header=TRUE, stringsAsFactors = F)

    # Validate metadata file
    miairr_metadata <- read.csv(miairr, sep="\t", stringsAsFactors = F)
    mandatory <- miairr_metadata[['Mandatory.BioSample.attribute']] == TRUE
    mandatory_fields <- setdiff(c(miairr_metadata[['AIRR.Formats.WG.field.name']][mandatory], "pcr_target_locus", "single_cell", "species"), "organism")

    if (!all(mandatory_fields %in% colnames(input))) {
        missing_fields <- mandatory_fields[mandatory_fields %in% colnames(input) == FALSE]
        if ("species" %in% missing_fields) {
            if ("organism" %in% colnames(input)) {
                warning("Missing 'species' field, using 'organism' field instead.")
                input[['species']] <- input[['organism']]
                missing_fields <- setdiff(missing_fields, "species")
            } else {
                stop("Missing mandatory field 'species'")
            }
        }
        warning("Missing MiAIRR fields: ", paste(missing_fields,collapse=", "))
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

    # Validate `species`
    species <- get_valid_species(valid_species, input[['species']])
    input[,'species'] <- species
    not_valid_species <- species[is.na(species)]

    input[,'valid_species'] <- TRUE
    if (length(not_valid_species)>0) {
        input[is.na(species),'valid_species'] <- FALSE
    }

    ## TODO: make a function to do column name validations
    # Validate collapseby
    collapseby <- strsplit(collapseby,",")[[1]]

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
    cloneby <- strsplit(cloneby,",")[[1]]
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


    # If reassign is false, check that tsv files have
    # all fields required by Immcantation (TODO)
    # and sequence_alignment has IMGT gaps
    # If reassign is true, sequence_alignment will be created later
    # in the pipeline, and will have imgt gaps
    if (reassign %in% c(TRUE, "true", "TRUE", T, "t", "T", 1, "1")) {
        input$valid_sequence_alignment <- TRUE
    } else {
        input$valid_sequence_alignment <- apply(input, 1, function(x) {
            if (x[['filetype']] == "tsv") {
                # if tsv files, look for presence of imgt gaps (.) in the
                # first rows
                any(grepl("\\.",read_rearrangement(x[['filename']], n_max=200)[['sequence_alignment']]))
            } else {
                # for non tsv files
                TRUE
            }
        })
    }

    check_fields <- c("valid_filename",
                        "valid_species",
                        "valid_collapseby",
                        "valid_cloneby",
                        "valid_single_cell",
                        "valid_pcr_target_locus",
                        "valid_sequence_alignment")
    not_valid_rows <- which(rowSums(input[, check_fields, drop=FALSE]) < length(check_fields))

    if (length(not_valid_rows) > 0) {
        message("Please review your input data. Invalid information found.")
    }

    list("table"=input,
        "validation_pass"=length(not_valid_rows) < 1
        )
}



# Settings
# Valid species. Use lower case.
valid_species <- list("human"=c("human", "homo sapiens", "h. sapiens", "hs"),
                        "mouse"=c("mouse", "mus musculus", "mm"))

# `species` must exist, then translate to name format
# suitable for igblast and downstream analysis
get_valid_species <- function(dictionary=valid_species, words) {
    dictionary <-  bind_rows(lapply(dictionary, function(x) {data.frame("species"=x)} ), .id="valid_species")
    sapply(words, function(word){
        if (is.na(word)) {
            warning("No `species` mapping available for: `",word,"`.")
            return (NA)
        }
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
