# Function to create the REcoTox project

#' @title Create project
#'
#' @description
#' This function initializes the \code{REcoTox} project in
#' the pre-defined \code{project_folder} in \code{project_path}.
#' It reads also the \code{EcoTox Knowledgebase} \emph{ascii}-files
#' provided in \code{database_path} and prepares the
#' \code{REcoTox} project.
#'
#' The prepared project is exported to \code{database_path} as default
#' project for re-use in other projects based on the same \code{EcoTox}
#' database version.
#'
#' @param database_path Path to the folder containing the expanded ascii files
#' of the current \code{EcoTox Knowledgebase}.
#'
#' @param project_path Path to the folder to store the processed project files.
#'
#' @param initalise_database_project Logical value to initialize a new default
#' \code{REcoTox} project in the \code{database_path}. The default value is
#' \code{FALSE} (values: \code{c(TRUE, FALSE)}).
#'
# @param initialise_project Logical value zu initialize the project folder in
# the \code{project_path}. The default value is \code{FALSE}
# (values:\code\{c(TRUE, FALSE)}).
#
#' @param load_default Logical value to load the default project in the
#' \code{database_path}. The default value is \code{FALSE} (values:
#' \code{c(TRUE, FALSE)}. The parameter is ignored, if the \code{initialise_
#' project} parameter is \code{TRUE}.
#'
#'
#' @author Tobias Schulze
#' @examples
#' # Initialize a new default \code{REcoTox} project in
#' # the \code{database_path} by reading the \code{
#' # EcoTox Knowledgebase} \emph{asii}-files.
#'
#' \dontrun{create_project(database_path = database_path,
#' project_path = project_path,
#' initalise_database_project = TRUE)}
#
#' @examples
#' # Initialize the \code{REcoTox} project folder in
#' # the \code{project_path} and read the default
#' # database from \code{database_path}.
#
#' \dontrun{create_project(database_path = database_path,
#' project_path = project_path,
#' initalise_project = TRUE,
#' load_default = TRUE)}
#
#' @export
#'

create_project <- function(database_path, project_path, initalise_database_project = FALSE,
                           initalise_project = FALSE, load_default = FALSE) {

  # initialise a new basic project, otherwise just load and update the project file
  if (initalise_database_project == TRUE) {
    project <- list()
    project$object <- list()
    project$database_path <- suppressWarnings(normalizePath(database_path))
    object <- project$object
    project$files <- normalizePath(list.files(database_path, pattern = "*.txt", full.names = T, recursive = T))
    message("[EcoToxR]:  Initializing a new project in the database folder. The old will be overwritten.")

    # read the file list and harmonise use of NA
    message("[EcoToxR]:  Read the Ecotox Knowledgebase ASCII files.")
    suppressWarnings({
      message("[EcoToxR]:  Read tests")
      object$tests <- readr::read_delim(file = project$files[grep("tests.txt", project$files)],
                                 delim = "|",
                                 na = c("NA","","NR","--","NC","/"),
                                 quote = "\"",
                                 show_col_types = FALSE)

      message("[EcoToxR]:  Read chemicals.")
      object$chemicals <- readr::read_delim(file = project$files[grep("chemicals.txt", project$files)],
                                     delim = "|",
                                     na = c("NA","","NR","--","NC","/"),
                                     quote = "\"",
                                     show_col_types = FALSE)

      object$chemprop <- tibble::tibble(create_chemical_properties(project$database_path))

      message("[EcoToxR]:  Read species.")
      object$species <- readr::read_delim(file = project$files[grep("species.txt", project$files)],
                                   delim = "|",
                                   na = c("NA","","NR","--","NC","/"),
                                   quote = "\"",
                                   show_col_types = FALSE)

      message("[EcoToxR]:  Read results.")

      object$results <- readr::read_delim(file = project$files[grep("results.txt", project$files)],
                                   delim = "|",
                                   na = c("NA","","NR","--","NC","/"),
                                   quote = "\"",
                                   show_col_types = FALSE)

      message("[EcoToxR]:  Read references.")
      object$references <- readr::read_delim(file = project$files[grep("references.txt", project$files)],
                                      delim = "|",
                                      na = c("NA","","NR","--","NC","/"),
                                      quote = "\"",
                                      show_col_types = FALSE)
    })

    #  Trim lists
    object$chemicals <- object$chemicals %>%
      dplyr::rename(compound_class = ecotox_group,
             dtxsid_ecotox = dtxsid)

    # harmonize name for query with chemicals
    #
    object$tests <- object$tests %>% dplyr::rename(cas_number = test_cas)

    # Add a column with a CAS in ususal format
    object$chemicals <- object$chemicals %>%
      dplyr::mutate(cas =
                  paste0(substr(as.character(object$chemicals$cas_number),
                  1,
                  nchar(as.character(object$chemicals$cas_number)) - 3),
                  "-",
                  substr(as.character(object$chemicals$cas_number),
                         nchar(as.character(object$chemicals$cas_number)) - 2,
                         nchar(as.character(object$chemicals$cas_number)) - 1),
                  "-",
                  substr(as.character(object$chemicals$cas_number),
                         nchar(as.character(object$chemicals$cas_number)),
                         nchar(as.character(object$chemicals$cas_number)))))

    message("[EcoToxR]:  Saving the default project in the database folder.")
    project$object <- object
    object <- NULL
    save(project,
         file = file.path(database_path, "project.RData"),
         compress = TRUE)
  }


  if (initalise_project == TRUE) {
    # Check the directory
    if (!dir.exists(project_path)) {
      message("[EcoToxR]:  Creating the project directory.")
      dir.create(project_path)
    } else if (!purrr::is_empty(project_path)) {
      message("[EcoToxR]:  The project directory is not empty. Please checkout content")
      message("[EcoToxR]:  and delete it or use another path before you continue.")
    } else if (purrr::is_empty(project_path)) {
      message("[EcoToxR]:  The project directory is empty.")
    }
    if (!exists("project") & load_default == FALSE) {
      load_default = TRUE
    }
    if (exists("project") & load_default == FALSE) {
      project$object$chemprop <-
        tibble::tibble(create_chemical_properties(project$database_path))
      project$project_path <- suppressWarnings(normalizePath(project_path))
    }
  }

  if (load_default == TRUE) {
    message("[EcoToxR]:  Loading the default project file.")
    load(file.path(database_path,"project.RData"))
    project$database_path <- suppressWarnings(normalizePath(database_path))
    project$project_path <- suppressWarnings(normalizePath(project_path))
    project$object$chemprop <- create_chemical_properties(database_path)
    project$files <-
      suppressWarnings(normalizePath(list.files(database_path,
                                                pattern = "*.txt",
                                                full.names = TRUE,
                                                recursive = TRUE)))
  }

  return(project)
}
