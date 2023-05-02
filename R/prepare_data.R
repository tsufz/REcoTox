# Prepares the data after import for downstream processing

#' @title Prepare data
#'
#' @description
#' This function prepares the data \code{REcoTox}. It can also handle
#' a previous created initial project, which can be stored in the end
#' of the preparation process to avoid re-processing data.
#'
#' If a project from a different folder is used, the current project folder
#' can be assigned.
#'
#' The \code{tests}, \code{chemicals}, \code{results}, \code{species}
#' and \code{references} imported from \code{ascii} files is merged.
#'
#' The function cleans the values from tilde, slash and asterics artefacts
#' and assign appropriate field types.
#'
#' The prepared project can be stored as initial project.
#'
#' @param project Name of the initial project in the \code{environment}.
#' The default value is \code{project}.
#'
#' @param load_initial_project Logical value to load an existing initial
#' project from the \code{project_folder} containing the pre-processed project.
#' The default value is \code{FALSE} (values: \code{c(TRUE, FALSE)}).
#'
#' @param new_project_path Path of the \code{project_folder}.
#' The default value is \code{NA}.
#'
#' @param reread_chemical_list Logical value if the chemical list (which
#' may be updated) is loaded in the re-loaded initial project.
#' The default value is \code{FALSE} (values:\code{c(TRUE, FALSE)}).
#
#' @param save_project Stores the initial project in the \code{project_folder}.
#' The default value is \code{TRUE} (values:\code{c(TRUE, FALSE)}).
#'
#' @author Tobias Schulze
#'
#' @examples
#' # Load the default \code{REcoTox} project in the \code{environment}
#' # and prepare the data for downstream processing.
#'
#' \dontrun{prepare_data(project = project,
#' load_initial_project = FALSE,
#' new_project_path = NA,
#' reread_chemical_list = FALSE,
#' save_project = FALSE)}
#
#' @examples
#'
#' # Load the default \code{REcoTox} project in the \code{environment}
#' # and prepare the data for downstream processing. Store the initial
#' # project in the \code{project_folder}.
#'
#' \dontrun{prepare_data(project = project,
#' load_initial_project = FALSE,
#' new_project_path = NA,
#' reread_chemical_list = FALSE,
#' save_project = TRUE)}
#'
#' @examples
#'
#' # Load the initial \code{REcoTox} project stored in the \code{project_folder}
#' # the \code{environment}, replace the \code{project_path}, add an updated
#' # \code{chemcial list} and the updated project in the \code{project_folder}.
#'
#' \dontrun{prepare_data(project = project,
#' load_initial_project = TRUE,
#' new_project_path = project_path,
#' reread_chemical_list = TRUE,
#' save_project = TRUE)}
#'
#' @export
#'
prepare_data <- function(project = project,
                         load_initial_project = FALSE,
                         new_project_path = NA,
                         reread_chemical_list = FALSE,
                         save_project = TRUE
                         ) {

  if (load_initial_project == TRUE) {
    initial_project <- normalizePath(file.path(project_path, "initial_project.RData"))
    if (file.exists(initial_project)) {
      message("[EcoToxR]:  Loading the initial project.")
      .tempenv <- new.env()
      load(file = initial_project, envir = .tempenv)
      if (!is.na(new_project_path)) {
        .tempenv$project$project_path <- normalizePath(project_path)
      }

      if (isTRUE(reread_chemical_list) ) {
          .tempenv$project <- readr::read_csv(file.path(database_path, "chemical_properties.csv"))
      }

      return(.tempenv$project)
      } else {
      message("[EcoToxR]:  The initial project does not exist. Please run 'Prepare_data' to prepare.")
      return()
    }
  }

  # Copy the object
  object <- project$object

  # record the settings
  object$parameters$save_project  <- save_project

  message("[EcoToxR]:  Merging tests, results and chemicals.")

  # Fix field type issues
  object$tests <- object$tests %>%
    dplyr::mutate(cas_number = as.double(cas_number),
                  reference_number = as.double(reference_number),
                  test_id = as.double(test_id)
                  )

  object$results <- object$results %>%
    dplyr::mutate(test_id = as.double(test_id),
                  result_id = as.double(result_id)
                  )
  object$references <- object$references %>%
    dplyr::mutate(reference_number = as.double(reference_number)
                  )

  # Merge the raw tables

  object$merged_results <- object$tests %>%
    dplyr::left_join(object$results, by = "test_id") %>%
    dplyr::left_join(object$chemicals, by = "cas_number") %>%
    dplyr::left_join(object$species, by = "species_number") %>%
    dplyr::left_join(object$references, by = "reference_number")

  # clean up data
  message("[EcoToxR]:  Cleaning up concentration, effects, endpoints, and measurements of asterics, tilde and slash artefacts. Fixing field type issues.")

  object$merged_results <- object$merged_results %>%
    dplyr::mutate(
      endpoint = as.character(gsub("\\*", "", endpoint)),
      endpoint = as.character(gsub("\\/", "", endpoint)),
      endpoint = as.character(gsub("\\~", "", endpoint)),

      measurement = as.character(gsub("\\*", "", measurement)),
      measurement = as.character(gsub("\\/", "", measurement)),
      measurement = as.character(gsub("\\~", "", measurement)),

      effect = as.character(gsub("\\*", "", effect)),
      effect = as.character(gsub("\\/", "", effect)),
      effect = as.character(gsub("\\~", "", effect)),

      # in some cases, the values are NA and thus trigger a warning
      conc1_mean = suppressWarnings(as.numeric(gsub("\\*", "", conc1_mean))),
      conc1_min = suppressWarnings(as.numeric(gsub("\\*", "", conc1_min))),
      conc1_max = suppressWarnings(as.numeric(gsub("\\*", "", conc1_max))),

      publication_year = as.integer(publication_year)
  )

  # Remove unnecessary tables to save space

  object$references <- NULL
  object$results <- NULL
  object$tests <- NULL

  # Do some last things and return
  project$object <- object

  if (save_project == TRUE) {
    message("[EcoToxR]:  Saving the initial project to the project folder.")
    file_name <- suppressWarnings(normalizePath(file.path(project$project_path, "initial_project.RData")))
    save(project, file = file_name, compress = TRUE)
  } else {

      return(project)
  }
  return(project)
}

