# Processing workflow to filter, amend and creation of the final longer
# data table.
#
#' @title Process data
#'
#' @description
#' This function processes the data in the \code{REcoTox} project according
#' to the filtering settings and creates a final longer data table.
#'
#' In brief, the \code{process_data} workflwo has four steps.
#'
#' Step 1: Read the filtering parameters, filters the data and exports
#' files for customized selection of endpoints and species.
#'
#' \itemize{
#' \item Read the parameters and store it in the project.
#' \item Filter the data according to given parameters.
#' \item Filter the data to volume related dosings.
#' \item Prepare and export the endpoint list for review.
#' \item Prepare and export the species list for review.
#' }
#'
#' Step 2: Read the reviewed endpoint and species tables, filter the data,
#' convert units to \code{mg/L} and prepare and export a list of compounds
#' with \code{mol/L} units for review and update of the \code{average_mass}
#' for unit transformation.
#'
#' \itemize{
#' \item Read the endpoint and species tables.
#' \item Filter the data according to selected endpoints and species.
#' \item Convert units to \code{mg/L}.
#' \item Prepare and export of \code{mol/L} units for review and amendment
#' of missing \code{average_mass} values.
#' }
#'
#' Step 3: Read the reviewed mol table and convert remaining \code{mol/L} units
#' to \code{mg/L}.
#' Prepare and export a list of compounds for review and amendment of
#' missing data.
#'
#' \itemize{
#' \item Read the revised mol table.
#' \item Convert remaining \code{mol/L} units to \code{mg/L}.
#' \item Prepare and export of the chemicals table for review.
#' }
#'
#' Step 4: Read the reviewed chemical table, filter excluded compounds,
#' update the basic chemical list in \code{database_folder}, build the final
#' longer list, convert \code{time} units to \code{h}, and export the list to the
#' \code{project_folder}.
#'
#' \itemize{
#' \item Read the revised chemical table.
#' \item Apply changes in the table.
#' \item Update the basic chemical list in the \code{database_folder} (optional).
#' \item Build the final results list in \code{longer} format.
#' \item Convert \code{day} units to \code{h}.
#' }
#'
#' @param project Name of the initial project in the \code{environment}.
#' The default value is \code{project}.
#'
#' @param dosing_group Dosings to be included in the analysis
#' Only the \code{water concentration} group is implemented.
#' The default value is \code{water_concentration}.
#'
#' @param duration_d Duration of the test in \code{days}.
#' The default value is \code{NA}. Typical values: \code{c("d", "dph", "dpf")}.
#'
#' @param duration_h Duration of the test in \code{hours}.
#' The default value is \code{NA}. Typical values: \code{c("h", "ht", "hph",
#' "hpf", "hbf", "hv")}.
#'
#' @param duration_m Duration of the test in \code{minutes}.
#' The default value is \code{NA}. Typical value: \code{mi}.
#'
#' @param ecotoxgroup Ecotoxgroup (species_group) to be selected for the
#' analysis.
#' The default value is \code{NA}. Typical values: \code{c(algae, crustaceans,
#' fish)}.
#'
#' @param effects Effects to included in the analysis.
#' For details, review the [US EPA Knowledgebase Code Appendix](https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix).
#' The default value is \code{NA}. Typical values: \code{c(algae, crustaceans,
#' fish)}.
#'
#' @param habitat Habitat to included in the analysis.
#' For details, review the [US EPA Knowledgebase Code Appendix](https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix).
#' The default value is \code{NA}. Allowed values: \code{c("Non-Soil","Water","Soil")}.
#'
##' @param kingdoms Kingdoms of algae to included in the analysis.
#' For details, review the [US EPA Knowledgebase Code Appendix](https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix).
#' The default value is \code{NA}.
#'
#' @param measurements Measurements to included in the analysis.
#' For details, review the [US EPA Knowledgebase Code Appendix](https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix).
#' The default value is \code{NA}.
#'
#' @param max_d Maximal duration of the test in \code{days}.
#' The default value is \code{NA}, otherwise the expected value is \code{d}.
#'
#' @param min_d Minimal duration of the test in \code{days}.
#' The default value is \code{NA} otherwise the expected value is \code{d}.
#'
#' @param max_h Maximal duration of the test in \code{hours}.
#' The default value is \code{NA} otherwise the expected value is \code{h}.
#'
#' @param min_h Minimal duration of the test in \code{hours}.
#' The default value is \code{NA}otherwise the expected value is \code{h}.
#'
#' @param max_m Maximal duration of the test in \code{minutes}.
#' The default value is \code{NA} otherwise the expected value is \code{m}.
#'
#' @param min_m Minimal duration of the test in \code{minutes}.
#' The default value is \code{NA} otherwise the expected value is \code{m}.
#'
#' @param remove_formulation Removes records related to formulations.
#' The default value is \code{NA}.
#'
#' @param sample_size Size of a sub-sample of records (for testing only).
#'
#' @param save_project_steps Stores the current project state in the \code{project_folder}.
#' The default value is \code{FALSE} (values: \code{c(FALSE, TRUE)}).
#'
#' @param species_selection Species to be selected for the analysis.
#' The default values is \code{NA}. Fore example: \code{c("all", "manual", "standard_test_species")}.
#'
#' @param state Set the processing state (for testing or re-run only).
#' The default value is \code{NA}.
#'
#' @param update_chemicals Update the basic chemical list in the
#' \code{database_folder}.
#' The default value is \code{TRUE} (values: \code{c(FALSE, TRUE)}).
#'
#' @author Tobias Schulze
#' @examples
#' # Run the workflow for processing step 1.
#'
#' \dontrun{project <- process_data(project = project,
#' duration_d = c("d", "dph", "dpf"),
#' duration_h = c("h", "ht", "hph", "hpf", "hbf", "hv"),
#' duration_m = mi,
#' ecotoxgroup = ecotoxgroup,
#' effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV"),
#' habitat = "Water",
#' max_d = 5,
#' min_d = 0,
#' max_h = 120,
#' min_h = 0,
#' max_m = 7200,
#' min_m = 0,
#' remove_formulation = FALSE,
#' save_project_steps = FALSE,
#' species_selection = "all")}
#
#' @examples
#' # Run the workflow for processing step 2.
#'
#' \dontrun{project <- process_data(project, save_project_steps = FALSE)}
#'
#' @examples
#' # Run the workflow for processing step 3.
#'
#' \dontrun{project <- process_data(project, save_project_steps = FALSE)}
#'
#' @examples
#' # Run the workflow for procssing step 4.
#'
#' \dontrun{project <- process_data(project, save_project_steps = FALSE, update_chemicals = FALSE)}
#'
#' @examples
#' # Run the workflow for procssing step 5.
#'
#' \dontrun{project <- aggregate_results(project = project, quantile = 0.05)}
#'
#' @export
#'
# functions

  process_data <- function(project = project,
                           dosing_group = "water_concentration",
                           duration_d = NA,
                           duration_h = NA,
                           duration_m = NA,
                           ecotoxgroup = NA,
                           effects = NA,
                           habitat = NA,
                           kingdoms = NA,
                           measurements = NA,
                           max_d = NA,
                           min_d = NA,
                           max_h = NA,
                           min_h = NA,
                           min_m = NA,
                           max_m = NA,
                           remove_formulation = FALSE,
                           sample_size = NA, # required for testing only
                           save_project_steps = FALSE,
                           species_selection = NA,
                           state = NA,
                           update_chemicals = FALSE
                           ){

    # Load data
    object <- project$object

    # Record the settings during first call (kept separated from next sections to maintain settings in this place)

    if (is.null(object$state) == TRUE){
        object$parameters$dosing_group <- dosing_group
        object$parameters$duration_h <- duration_h
        object$parameters$duration_d <- duration_d
        object$parameters$duration_m <- duration_m
        object$parameters$effects <- effects
        object$parameters$ecotoxgroup <- ecotoxgroup
        object$parameters$habitat <- habitat
        object$parameters$measurements <- measurements
        object$parameters$kingdoms <- kingdoms
        object$parameters$max_d <- max_d
        object$parameters$min_d <- min_d
        object$parameters$max_h <- max_h
        object$parameters$min_h <- min_h
        object$parameters$max_m <- max_m
        object$parameters$min_m <- min_m
        object$parameters$remove_formulation <- remove_formulation
        object$parameters$sample_size <- sample_size
        object$parameters$species_selection <- species_selection
        object$parameters$update_chemicals <- update_chemicals
    }

    # Workflow step 1
    if (is.null(object$state) == TRUE) {

      # select only water organisms and merge tables
      message("[EcoToxR]:  Subsetting to habitat.")
      object$results_filtered <- object$merged_results %>%
        dplyr::filter(organism_habitat %in% habitat)

      # Do preliminary filtering (group and effects)
      if (!is.na(ecotoxgroup)) {
          message("[EcoToxR]:  Filtering the datasets by the ecotox group (species group).")
          object$results_filtered <-
            object$results_filtered %>%
            dplyr::filter(ecotox_group %like% ecotoxgroup)
      }

      # Too complicated, needs refactorization see issue #9 https://git.ufz.de/WANA/REcoTox/-/issues/9
      if (is.na(dosing_group)) {
          message("[EcoToxR]:  The dosing group cannot be empty or set to \"NA\"")
      } else if (dosing_group == "water_concentration"){
          message("[EcoToxR]:  Reducing to mass/water related units.")
          unit_pattern <- "(^(|(A|a)(e|i|E|I) )(|m|n|u|p)(g|mol)/(|m|u|p|d|)(L|l|m3)$)|(^(p)(p)(m|t|b)$)"
          object$results_filtered <- object$results_filtered %>%
            dplyr::filter(conc1_unit %like% unit_pattern)

          message("[EcoToxR]:  The following units have been found in the data:")
          print(dplyr::pull(object$results_filtered %>%
                       dplyr::select(conc1_unit) %>%
                       unique() %>%
                       dplyr::group_by(conc1_unit) %>%
                       dplyr::arrange(.by_group = TRUE)))
      } else if (!dosing_group == "water_concentration") {
          message("[EcoToxR]:  The dosing group must be set to \"water_concentration\"")
      }


      # Subset algae data to kingdoms
      if (ecotoxgroup == "Algae" && !is.na(kingdoms)) {
          object$results_filtered <- object$results_filtered %>%
            dplyr::filter(kingdom %in% kingdoms)
      }


      # Select by effect
      if (!is.logical(effects)) {
          message("[EcoToxR]:  Filtering the datasets by the effects")
          object$results_filtered <- object$results_filtered %>%
            dplyr::filter(effect %in% effects)

          message("[EcoToxR]:  The following effects are available:")
          print(pull(object$results_filtered %>%
                       dplyr::select(effect) %>%
                       unique() %>%
                       dplyr::group_by(effect) %>%
                       dplyr::arrange(.by_group = TRUE)))
      }

      # Select by measurement
      if (!is.logical(measurements)) {
          message("[EcoToxR]:  Filtering the datasets by the measurement.")
          object$results_filtered <- object$results_filtered %>%
            dplyr::filter(measurement %in% measurements)

          message("[EcoToxR]:  The following measurements are available:")
          print(pull(object$results_filtered %>%
                       dplyr::select(measurement) %>%
                       unique() %>%
                       dplyr::group_by(measurement) %>%
                       dplyr::arrange(.by_group = TRUE)))
      }

      # Remove formulations
      if (remove_formulation == TRUE) {
          message("[EcoToxR]:  Removing formulations.")
          object$results_filtered <- object$results_filtered %>%
            dplyr::filter(conc1_type == "A")
          removed_formulation <- object$results_filtered %>%
            dplyr::filter(conc1_type == "T")
          message(paste0("[EcoToxR]:  The dataset includes ",
                         nrow(object$results_filtered),
                         " records of active ingredients."))

          message(paste0("[EcoToxR]:  ",
                         nrow(removed_formulation),
                         " record(s) of formulations was/were removed."))

          rm(removed_formulation)
      }



    # Handling durations
      # Filter the time groups (minutes, hours, days)
    results_per_days <- object$results_filtered %>%
      dplyr::filter(obs_duration_unit %in% duration_d)

    results_per_hours <- object$results_filtered %>%
      dplyr::filter(obs_duration_unit %in% duration_h)

    results_per_minutes <- object$results_filtered %>%
      dplyr::filter(obs_duration_unit %in% duration_m)


    # Subset to minimal and maximal durations

    results_per_days <- results_per_days %>%
      dplyr::filter(obs_duration_mean >= min_d)

    results_per_days <- results_per_days %>%
      dplyr::filter(obs_duration_mean <= max_d)

    results_per_hours <- results_per_hours %>%
      dplyr::filter(obs_duration_mean >= min_h)

    results_per_hours <- results_per_hours %>%
      dplyr::filter(obs_duration_mean <= max_h)

    results_per_minutes <- results_per_minutes %>%
      dplyr::filter(obs_duration_mean >= min_m)

    results_per_minutes <- results_per_minutes %>%
      dplyr::filter(obs_duration_mean <= max_m)

    object$results_filtered <-
      dplyr::bind_rows(results_per_days,
                       results_per_hours,
                       results_per_minutes)

    rm(results_per_days, results_per_hours, results_per_minutes)


    # Prepare and export a list of species
    # species_selection = c("all", "selected", "standard_test_species")
    if (is.na(species_selection) | species_selection == "all") {
        species <- object$results_filtered %>%
          dplyr::select(species_number, ecotox_group) %>%
          unique() %>%
          dplyr::select(species_number) %>%
          dplyr::left_join(object$species, by = "species_number")

        species <- species %>%
          tibble::add_column("include_species" = 1, .before = 1) %>%
          tibble::add_column("count_of_records" = 0, .after = 1)

        for (i in 1:nrow(species)) {
            species$count_of_records[i] <-
              length(which(object$results_filtered[, "species_number"] ==
                                                            species[i, "species_number"][[1]]))
        }

        species <- species %>%
          dplyr::arrange(-count_of_records)

        message("[EcoToxR]:  The table with the species was written to the project folder.")
        #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
        readr::write_csv(species, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")))))
    } else if (species_selection == "manual") {

      species <- object$results_filtered %>%
        dplyr::select(species_number, ecotox_group) %>%
        unique() %>%
        dplyr::select(species_number) %>%
        dplyr::left_join(object$species, by = "species_number") %>%
        tibble::add_column("include_species" = 0, .before = 1) %>%
        tibble::add_column("count_of_records" = 0, .after = 1)


        for (i in 1:nrow(species)) {
          species$count_of_records[i] <-
            length(which(object$results_filtered[, "species_number"] ==
                                                                    species[i, "species_number"][[1]]))
        }

        species <- species %>% dplyr::arrange(-count_of_records)

        message("[EcoToxR]:  The table with the species was written to the project folder.")
        #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
       readr::write_csv(species, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      } else if (species_selection == "standard_test_species") {

          species <- object$results_filtered %>%
            dplyr::select(species_number, ecotox_group) %>%
            dplyr::filter(ecotox_group %like% "Standard Test Species") %>%
            unique() %>%
            dplyr::select(species_number) %>%
            dplyr::left_join(object$species, by = "species_number")

          species <- species %>%
            tibble::add_column("include_species" = 1, .before = 1) %>%
            tibble::add_column("count_of_records" = 0, .after = 1)

          for (i in 1:nrow(species)) {
              species$count_of_records[i] <-
                length(which(object$results_filtered[, "species_number"] ==
                                                              species[i, "species_number"][[1]]))
          }
          species <- species %>%
            dplyr::arrange(-count_of_records)

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          readr::write_csv(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      } else {
          species <- object$results_filtered %>%
            dplyr::select(species_number, ecotox_group) %>%
            unique() %>%
            dplyr::select(species_number) %>%
            dplyr::left_join(object$species, by = "species_number")

          species <- species %>%
            tibble::add_column("include_species" = 1, .before = 1) %>%
            tibble::add_column("count_of_records" = 0, .after = 1)

          for (i in 1:nrow(species)) {
              species$count_of_records[i] <-
                length(which(object$results_filtered[, "species_number"] ==
                                                              species[i, "species_number"][[1]]))
          }
          species <- species %>%
            dplyr::arrange(-count_of_records)

          message("[EcoToxR]:  Species selection has unvalid input. The table with the species was written to the project folder for manual review.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          readr::write_csv(species, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      }

    # Select ecotox/species group (e.g. "Algae")
    message("[EcoToxR]:  The following ecotox/species group(s) are included:")
    print(dplyr::pull(species %>%
                        dplyr::select(ecotox_group) %>%
                        unique() %>%
                        dplyr::group_by() %>%
                        dplyr::arrange(.by_group = TRUE)))


      endpoints <- object$results_filtered %>%
        dplyr::select(endpoint) %>%
        unique() %>%
        tibble::add_column("include_endpoint" = 1, .before = 1) %>%
        tibble::add_column("count_of_records" = 0, .after = 1)

      for (i in 1:nrow(endpoints)) {
          endpoints$count_of_records[i] <-
            length(which(object$results_filtered[, "endpoint"] ==
                                                            endpoints[i, "endpoint"][[1]]))
      }

      endpoints <- endpoints %>%
        dplyr::group_by(endpoint) %>%
        dplyr::arrange(.by_group = TRUE)

      message("[EcoToxR]:  The table with the endpoints was written to the project folder.")
      message("[EcoToxR]:  Edit the file(s) and re-run the workflow.")

      write_csv(endpoints, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_endpoint_selection.csv")))))

      object$state <- 1

      project$object <- object
      save_project(project, save_project_steps)

      #äfwrite(object$mortality,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_mortality_unfiltered.csv")))),  sep = ",", dec = ".")
      return(project)
    }

    # Workflow step 2

    if(object$state == 1){
    ecotoxgroup <- object$parameters$ecotoxgroup
    # Filter the selected endpoints
    endpoints <- readr::read_csv(file.path(project_path, paste0(tolower(ecotoxgroup), "_endpoint_selection.csv")), show_col_types = FALSE)

    if(!all(unique(endpoints$include_endpoint) %in% c(0, 1))) {
        stop("[EcoToxR]:  Only 0 or 1 is allowed in the endpoint selection filter. Check the endpoint selection list and re-run the workflow.")
    }

    object$results_filtered <- object$results_filtered %>%
      dplyr::left_join(endpoints %>%
                       dplyr::select(include_endpoint, endpoint),
                       by = "endpoint") %>%
      dplyr::filter(include_endpoint == 1)

    # Filter the selected species

    species <- read_csv(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")), show_col_types = FALSE)

    object$results_filtered <- object$results_filtered %>%
        left_join(species %>%
        select(include_species, species_number), by = "species_number") %>%
        filter(include_species == 1)

    # Convert the units
    object$results_filtered <- convert_units(object$results_filtered, sample_size)

    # extract mol related records and export for review, alternatively export the file for exclusion review
    if (nrow(object$results_filtered %>% select(conc1_unit) %>% filter(conc1_unit %like% "mol/L")) > 0){
        object <- export_mol_units(object = object)
        object$state <- 2
        project$object <- object
        save_project(project, save_project_steps)
    return(project)
    } else {

    #object <- export_exclude_list(object,project$project_path)
    export_chemical_list(object, project$project_path)
    object$state <- 3
    project$object <- object
    save_project(project, save_project_steps)
    return(project)
    }

}

    # workflow step 3

    if(object$state == 2){
      ecotoxgroup <- object$parameters$ecotoxgroup
      # update records with mol related concentrations
      if (file.exists(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mol_weight.csv")))){
        object <-
          update_mol_units(object = object,
                           database_path = project$database_path,
                           project_path = project$project_path)
        export_chemical_list(object, project$project_path)
        object$state <- 3
        project$object <- object
        save_project(project, save_project_steps)
        return(project)
      }
    }

    # Workflow step 4

    if(object$state == 3){

      # This routine just exports an exclusion and inclusion list, but does not remove chemicals from chemical list
      object <- remove_excluded_chemicals(object, project$project_path)

      if(update_chemicals == TRUE){
        object <-
          update_chemical_list(object,
                               project_path = project_path,
                               database_path = database_path)
      }
      object <- build_final_list(object)
      object <- calculate_hours(object)

      object$state <- 4
      readr::write_csv(x = object$results,file = suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_final_results.csv")))))
      project$object <- object
      save_project(project, save_project_steps)
      message("[EcoToxR]:  The data pre-processing is finalised.")
      message("[EcoToxR]:  Please run the pivot export workflow.")
      return(project)
    }
  }
