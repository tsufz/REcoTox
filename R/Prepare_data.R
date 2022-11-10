#' @export
#'


prepare_data <- function(project,
                         load_initial_project = FALSE,
                         new_project_path = NA,
                         reread_chemical_list = NA,
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

      if ( isTRUE(reread_chemical_list) ) {
          .tempenv$project <- read_csv(file.path(database_path, "chemical_properties.csv"))
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
  object$chemicals <- suppressWarnings(object$chemicals %>% mutate(cas_number = as.integer(cas_number)))
  object$tests <- suppressWarnings(object$tests %>% mutate(cas_number = as.integer(cas_number), reference_number = as.integer(reference_number), test_id = as.integer(test_id)))
  object$results <- suppressWarnings(object$results %>% mutate(test_id = as.integer(test_id)))
  object$references <- suppressWarnings(object$references %>% mutate(reference_number = as.integer(reference_number)))

  # Merge the raw tables

  object$merged_results <- object$tests %>% left_join(object$results, by = "test_id") %>%
      left_join(object$chemicals, by = "cas_number") %>%
      left_join(object$species, by = "species_number") %>%
      left_join(object$references, by = "reference_number")

  # clean up data
  message("[EcoToxR]:  Cleaning up concentration, effects, endpoints, and measurements of asterics, tilde and slash artefacts. Fixing field type issues.")

  object$merged_results <- object$merged_results %>%
      mutate(endpoint = as.character(gsub("\\*", "", endpoint)),
             endpoint = as.character(gsub("\\/", "", endpoint)),
             endpoint = as.character(gsub("\\~", "", endpoint)),

             measurement = as.character(gsub("\\*", "", measurement)),
             measurement = as.character(gsub("\\/", "", measurement)),
             measurement = as.character(gsub("\\~", "", measurement)),

             effect = as.character(gsub("\\*", "", effect)),
             effect = as.character(gsub("\\/", "", effect)),
             effect = as.character(gsub("\\~", "", effect)),

             conc1_mean = suppressWarnings(as.numeric(gsub("\\*", "", conc1_mean))),
             conc1_min = suppressWarnings(as.numeric(gsub("\\*", "", conc1_min))),
             conc1_max = suppressWarnings(as.numeric(gsub("\\*", "", conc1_max))),

             publication_year = suppressWarnings(as.integer(publication_year))
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

