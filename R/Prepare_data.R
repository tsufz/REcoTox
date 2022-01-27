prepare_data <- function(project,
                         effects = c("MOR","GRO","DEV"),
                         habitat = c("Non-Soil","Water","Soil"),
                         remove_formulation = FALSE,
                         save_project = TRUE,
                         new_project_path = NA,
                         load_initial_project = FALSE) {

  if (load_initial_project == TRUE) {
    initial_project <- normalizePath(file.path(project_path, "initial_project.RData"))
    if (file.exists(initial_project)) {
      message("[EcoToxR]:  Loading the initial project.")
      .tempenv <- new.env()
      load(file = initial_project, envir = .tempenv)
      if (!is.na(new_project_path)) {
        .tempenv$project$project_path <- normalizePath(project_path)
      }
      return(.tempenv$project)
      } else {
      message("[EcoToxR]:  The initial project does not exist. Please run 'Prepare_data' to prepare.")
      return()
    }
  }
  object <- project$object

  # record the settings
  object$parameters$habitat <- habitat
  object$parameters$effects <- effects
  object$parameters$remove_formulation <- remove_formulation
  object$parameters$save_project  <- save_project

  # select only water organisms and merge tables
  message("[EcoToxR]:  Subsetting to habitat.")
  object$tests_habitat <- object$tests %>% filter(organism_habitat %in% habitat)

  # Merge the raw tables
  message("[EcoToxR]:  Merging tests, results and chemicals.")
  object$tests_results <- object$tests_habitat %>% left_join(object$results, by = "test_id")
  object$tests_results <- object$tests_results %>% mutate(cas_number = as.integer(cas_number), endpoint = as.character(endpoint))
  object$chemicals <- suppressWarnings(object$chemicals %>% mutate(cas_number = as.integer(cas_number)))

  object$tests_results_chemicals <- tibble(object$tests_results %>% left_join(object$chemicals, by = "cas_number"))
  object$test_results_all <- object$tests_results_chemicals %>% left_join(object$species, by = "species_number")
  object$test_results_all <- object$test_results_all %>% left_join(object$references, by = "reference_number")

  # delete intermediate tables
  object$tests <- NULL
  #object$chemicals <- NULL
  #object$species <- NULL
  object$results <- NULL
  object$tests_habitat <- NULL
  object$tests_results <- NULL
  object$tests_results_chemicals <- NULL
  object$references <- NULL

  # clean up data
  message("[EcoToxR]:  Cleaning up effects, endpoints and measurements of asterics. tilde and slash artefacts.")

  object$test_results_all <- object$test_results_all %>% mutate(endpoint = as.character(gsub("[*].*$", "", endpoint)))
  object$test_results_all <- object$test_results_all %>% mutate(endpoint = as.character(gsub("[/].*$", "", endpoint)))
  object$test_results_all <- object$test_results_all %>% mutate(measurement = as.character(gsub("[*].*$", "", measurement)))
  object$test_results_all <- object$test_results_all %>% mutate(measurement = as.character(gsub("[/].*$", "", measurement)))
  object$test_results_all <- object$test_results_all %>% mutate(effect = as.character(gsub("[*].*$", "", effect)))
  object$test_results_all <- object$test_results_all %>% mutate(effect = as.character(gsub("[/].*$", "", effect)))

  # Select data by effects
  message(paste0("[EcoToxR]:  Shrinking to effects: ", paste(effects, collapse = ", "), "."))
  object$parameters$effects <- effects

  object$all_selected_effects <- object$test_results_all %>% filter(effect %in% effects)

  message("[EcoToxR]:  The following effects are included:")
  print(unique(object$all_selected_effects$effect))
  object$test_results_all <- NULL

  if (remove_formulation == TRUE) {
    message("[EcoToxR]:  Removing formulations.")
    object$all_selected_effects <- object$all_selected_effects[conc1_type %like% "A"]
    object$removed_formulation <- object$all_selected_effects[conc1_type %like% "T"]
    message(paste0("[EcoToxR]:  The dataset includes ", nrow(object$all_selected_effects), " records of active ingredients."))
    message(paste0("[EcoToxR]:  ", nrow(object$removed_formulation), " record(s) of formulations was/were removed."))
  }

  message("[EcoToxR]:  Reducing to mass/water related units.")
  unit_pattern <- "(^(|(A|a)(e|i|E|I) )(|m|n|u|p)(g|mol)/(|m|u|p|d|)(L|l|m3)$)|(^(p)(p)(m|t|b)$)"
  object$all_selected_effects <- object$all_selected_effects %>% filter(conc1_unit %like% unit_pattern)

  message("[EcoToxR]:  The following units have been found in the data:")
  print(unique(object$all_selected_effects$conc1_unit))

  object$all_selected_effects <-  object$all_selected_effects %>% drop_na(endpoint)

  object$all_selected_effects <- data.table(reduce(list(
  object$chemicals[, colnames(object$chemicals) %in% c("cas_number", "chemical_name", "dtxsid"), with = FALSE],
  object$all_selected_effects
  ),
  left_join,
  by = "cas_number"))

  colnames(object$all_selected_effects)[2] <- "chemical_name"

  object$tests_habitat <- NULL
  object$tests_results <- NULL
  object$test_results_all <- NULL
  project$object <- object

  if (save_project == TRUE) {
    message("[EcoToxR]:  Saving the initial project to the project folder.")
    file_name <- suppressWarnings(normalizePath(file.path(project$project_path, "initial_project.RData")))
    save(project, file = file_name, compress = TRUE)
  }
  return(project)
}
