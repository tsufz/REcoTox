## ----biocstyle, echo = FALSE, results = "asis"--------------------------------
BiocStyle::markdown()

## ----init, message = FALSE, echo = FALSE, results = "hide"--------------------
## Silently loading all packages
library(BiocStyle)
library(desc)
library(kableExtra)
library(tidyverse)

## ----load REcoTox package, eval = FALSE, echo = TRUE, message = FALSE, warning = FALSE----
#  # Load the REcoTox package
#  library(REcoTox)

## ----R Documentation, echo = TRUE, eval = FALSE-------------------------------
#  # Documentation of REcoTox
#  help(package = "REcoTox")

## ----initialize folders, echo = TRUE, message = FALSE, warning = FALSE, eval = TRUE----
# Path of the project folder
project_folder <- "REcoTox_demo"

database_folder <- system.file("extdata/database_folder", package="REcoTox")
# The project folder is created in the home directory
project_path <- normalizePath(ifelse(.Platform$OS.type == "unix",
    paste0("~/", project_folder),
    paste0(
        Sys.getenv("HOMEPATH"),
        "\\",
        project_folder
    )
))

# An existing folder is deleted
#if (dir.exists(project_folder)) {
#    unlink(project_folder, recursive = TRUE)
#}

## ----create project, echo = TRUE, message = FALSE, warning = FALSE, eval = TRUE----
project <- REcoTox::create_project(database_path = database_folder,
                          project_path,
                          initalise_database_project = TRUE, # create the basic project from current ASCII files in DB folder
                          initalise_project = TRUE, # initializes the project folder
                          load_default = FALSE) # loads the default project in the project folder in the memoryfault_example = TRUE

file.copy(
                    from = system.file(
                        "extdata",
                        "Query_EcoTox_DB.R",
                        package = "REcoTox"
                    ),
                    to = normalizePath(
                        path = file.path(
                            project_folder,
                            "Query_EcoTox_DB.R"
                        ),
                        winslash = "\\",
                        mustWork = FALSE
                    ),
                    overwrite = TRUE
                )


## ----list project folder------------------------------------------------------
# List files and directories in project_folder
list.files(project_folder, recursive = TRUE, include.dirs = TRUE)

## ----list database folder-----------------------------------------------------
# List files and directories in project_folder
list.files(database_folder, recursive = TRUE, include.dirs = TRUE)

## ----view chemical_properties, echo = TRUE, eval = TRUE, message = TRUE-------
# Review of the chemical properties
chemical_properties <- readr::read_csv(file = normalizePath(path = file.path(
    database_folder,
    "chemical_properties.csv"
), ), show_col_types = FALSE)

kable(
    chemical_properties %>%
        select(cas_number:dtxsid_ecotox) %>% 
        head(5),
    format = "latex", digits = 2
)

## ----view results, echo = TRUE, eval = TRUE, message = TRUE-------------------
# Review of the result table
results <-
    readr::read_delim(
        file = normalizePath(
            path = file.path(
                database_folder,
                "results.txt"
            ),
        ),
        show_col_types = FALSE,
        delim = "|"
        
    )

kable(
    results  %>%
        select(result_id:sample_size_mean) %>% 
        head(5),
    format = "latex", digits = 2
)

## ----view chemicals, echo = TRUE, eval = TRUE, message = TRUE-----------------
# Review of the substance_table
substances <-
    readr::read_delim(
        file = normalizePath(
            path = file.path(
                database_folder,
                "validation",
                "chemicals.txt"
            ),
        ),
        show_col_types = FALSE,
        delim = "|"
        
    )

kable(
    substances %>%
        select(cas_number:ecotox_group) %>% 
        head(5),
    format = "latex", digits = 2
)

## ----view references, echo = TRUE, eval = TRUE, message = TRUE----------------
# Review of the substance_table
references <-
    readr::read_delim(
        file = normalizePath(
            path = file.path(
                database_folder,
                "validation",
                "references.txt"
            ),
        ),
        show_col_types = FALSE,
        delim = "|"
        
    )

kable(
    references %>%
        select(reference_number:author) %>% 
        head(5),
    format = "latex", digits = 2
)

## ----view species, echo = TRUE, eval = TRUE, message = TRUE-------------------
# Review of the substance_table
species <-
    readr::read_delim(
        file = normalizePath(
            path = file.path(
                database_folder,
                "validation",
                "species.txt"
            ),
        ),
        show_col_types = FALSE,
        delim = "|"
        
    )

kable(
    species %>%
        select(species_number:kingdom) %>% 
        head(5),
    format = "latex", digits = 2
)

## ----initialize databases, echo = TRUE, message = FALSE, warning = FALSE, eval = FALSE----
#  project <- REcoTox::create_project(database_path = database_folder,
#                            project_path,
#                            initalise_database_project = TRUE,
#                            initalise_project = TRUE,
#                            load_default = FALSE)
#  )

## ----initialize project, echo = TRUE, message = FALSE, warning = FALSE, eval = FALSE----
#  project <- REcoTox::prepare_data(project = project,
#                          load_initial_project = FALSE,
#                          new_project_path = NA,
#                          save_project = TRUE
#  )

## ----sessioninfo, echo = TRUE, eval = TRUE, message = FALSE-------------------
sessionInfo()

## ----clean_up, echo = FALSE, results = "asis", eval = FALSE-------------------
#  #unlink(project_folder, recursive = TRUE)

