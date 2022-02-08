# functions

  process_data <- function(project,
                           ecotoxgroup = NA,
                           sample_size = NA,
                           species_selection = "selected",
                           measurements = NA,
                           min_h = 0, min_d = 0,
                           max_h = 120, max_d = 5,
                           duration_h = c("h","ht","hph","hpf","hv","hbf"),
                           duration_d = c("d","dph","dpf"),
                           kingdoms = c("Chromista","Plantae","Monera"),
                           water_solubility_cut_off = FALSE,
                           water_solubility_threshold = NA,
                           remove_operator_records = FALSE,
                           save_project_steps = TRUE,
                           update_chemicals = FALSE,
                           all_species = FALSE
                           ){
    object <- project$object


    # Record the settings
    object$parameters$ecotoxgroup <- ecotoxgroup
    object$parameters$measurements <- measurements
    object$parameters$min_h <- min_d
    object$parameters$min_d <- min_d
    object$parameters$max_h <- max_d
    object$parameters$max_d <- max_h
    object$parameters$duration_h <- duration_h
    object$parameters$duration_d <- duration_d
    object$parameters$kingdom <- kingdoms
    object$parameters$sample_size <- sample_size
    object$parameters$species_selection <- species_selection
    object$parameters$all_species <- all_species

    # Workflow step 1
    if (is.null(object$state) == TRUE) {
      # Do preliminary filtering (group and effects)
      if (!is.na(ecotoxgroup)) {
          message("[EcoToxR]:  Filtering the datasets by the ecotox group (species group).")
          object$mortality <- dplyr::tibble(object$all_selected_effects %>% filter(ecotox_group %like% ecotoxgroup))
      } else {
          object$mortality <- dplyr::tibble(object$all_selected_effects)
       }

      object$all_selected_effects <- NULL

      # Select mortality records only
      message("[EcoToxR]:  The following ecotox/species group(s) are available:")
      print(unique(object$mortality$ecotox_group))

      if (!is.na(measurements)) {
        message("[EcoToxR]:  Filtering the datasets by the measurement.")
        object$mortality <- dplyr::tibble(object$mortality %>% dplyr::filter(measurement %in% measurements))

        message("[EcoToxR]:  The following measurements are available:")
        print(sort(unique(object$mortality$measurement)))
      }

      # Subset algae data to kingdom
      if (ecotoxgroup == "Algae" && !is.na(kingdoms)) {
        object$mortality <- object$mortality %>% dplyr::filter(kingdom %in% kingdoms)
      }


      # Select given durations
      mortality_h <- object$mortality %>% dplyr::filter(obs_duration_unit %in% duration_h)
      mortality_d <- object$mortality %>% dplyr::filter(obs_duration_unit %in% duration_d)

      # Subset to minimal and maximal durations

      mortality_h <- mortality_h %>% dplyr::filter(obs_duration_mean >= min_h)
      mortality_h <- mortality_h %>% dplyr::filter(obs_duration_mean <= max_h)

      mortality_d <- mortality_d %>% dplyr::filter(obs_duration_mean >= min_d)
      mortality_d <- mortality_d %>% dplyr::filter(obs_duration_mean <= max_d)

      object$mortality <- dplyr::tibble(rbind(mortality_h, mortality_d))

      mortality_h <- NULL
      mortality_d <- NULL

      # Prepare and export a list of endpoints
      endpoints <- dplyr::tibble(object$mortality %>%
                                     select(endpoint) %>%
                                     unique())

      # Prepare and export a list of species
      # species_selection = c("all", "selected", "standard_test_species")
      if (species_selection == "selected") {

          species <- object$mortality %>% dplyr::select(species_number, ecotox_group) %>%
              unique() %>%
              dplyr::select(species_number) %>%
              dplyr::left_join(object$species, by = "species_number") %>%
              tibble::add_column("include_species" = 0, .before = 1) %>%
              tibble::add_column("count_of_records" = 0, .after = 1)


          for (i in 1:nrow(species)) {
            species$count_of_records[i] <- length(which(object$mortality[, "species_number"] ==
                                                                    species[i, "species_number"][[1]]))
          }

          species <- species %>% arrange(-count_of_records)

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          fwrite(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup),"_species_selection.csv")))),  sep = ",", dec = ".")
      } else if (species_selection == "standard_species") {

          species <- object$mortality %>% select(species_number, ecotox_group) %>%
              filter(ecotox_group %like% "Standard Test Species") %>%
              unique() %>%
              select(species_number) %>%
              left_join(object$species, by = "species_number")

          species <- species %>%
              tibble::add_column("include_species" = 1, .before = 1) %>%
              tibble::add_column("count_of_records" = 0, .after = 1)

          for (i in 1:nrow(species)) {
              species$count_of_records[i] <- length(which(object$mortality[, "species_number"] ==
                                                              species[i, "species_number"][[1]]))
          }
          species <- species %>% arrange(-count_of_records)

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          fwrite(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))),  sep = ",", dec = ".")

      }

      endpoints <- object$mortality %>%
          select(endpoint) %>%
          unique() %>%
          tibble::add_column("include_endpoint" = 1, .before = 1) %>%
          tibble::add_column("count_of_records" = 0, .after = 1)

      for (i in 1:nrow(endpoints)) {
          endpoints$count_of_records[i] <- length(which(object$mortality[, "endpoint"] ==
                                                            endpoints[i, "endpoint"][[1]]))
      }

      endpoints <- endpoints %>% arrange(endpoint)

      message("[EcoToxR]:  The table with the endpoints was written to the project folder.")
      message("[EcoToxR]:  Edit the file(s) and re-run the workflow.")
      fwrite(endpoints, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_endpoint_selection.csv")))),  sep = ",", dec = ".")

      object$state <- 1

      project$object <- object
      save_project(project, project_path, save_project_steps)

      fwrite(object$mortality,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_mortality_unfiltered.csv")))),  sep = ",", dec = ".")
      return(project)
    }

    # Workflow step 2

    if(object$state == 1){
      endpoints <- dplyr::tibble(fread(file.path(project_path, paste0(tolower(ecotoxgroup), "_endpoint_selection.csv"))))
      if(!all(unique(endpoints$include_endpoint) %in% c(0, 1))) {
        stop("[EcoToxR]:  Only 0 or 1 is allowed in the endpoint selection filter. Check the endpoint selection list and re-run the workflow.")
      }

      if(unique(endpoints$include_endpoint) == 0) {
        stop("[EcoToxR]:  There are only 0 in the endpoint selection filter. Check the species selection list and re-run the workflow.")
      }

      object$mortality <- object$mortality %>% left_join(endpoints %>% select(include_endpoint, endpoint), by = "endpoint") # merge

      # Tidy concentration field
      #
      object$mortality <- tidy_conc_values(object$mortality)
      object$mortality_filtered <- object$mortality %>% filter(include_endpoint == 1)
      object$mortality_filtered_removed_endpoint <- object$mortality %>% filter(is.na(include_endpoint) || include_endpoint == 0)

      if(object$parameters$all_species == FALSE) {
        species <- dplyr::tibble(fread(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")), sep = ",", dec = "."))
        if(!all(unique(species$include_species) %in% c(0, 1))) {
          stop("[EcoToxR]:  Only 0 or 1 is allowed in the species selection filter. Check the species selection list and re-run the workflow.")
        }
        if(unique(species$include_species) == 0) {
          stop("[EcoToxR]:  There are only 0 in the species selection filter. Check the species selection list and re-run the workflow.")
        }

        object$mortality_filtered <- object$mortality_filtered %>%
            left_join(species %>%
            select(include_species, species_number), by = "species_number") # merge

        object$mortality_filtered <- object$mortality_filtered %>% filter(include_species == 1)
        object$mortality_filtered_removed_endpoint <- object$mortality %>%
            left_join(species %>% select(include_species, species_number), by = "species_number") %>%
            filter(is.na(include_species) || include_species == 0)


        fwrite(object$mortality_filtered, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mortality_filtered.csv")))),  sep = ",", dec = ".")
        fwrite(object$mortality_filtered, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mortality_removed_species.csv")))),  sep = ",", dec = ".")
      }

      # Check measurements

      # Convert the units and save the file
      object$mortality_filtered <- convert_units(object$mortality_filtered, sample_size)
      fwrite(object$mortality_filtered_removed_endpoint, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mortality_removed_endpoints.csv")))),  sep = ",", dec = ".")

      # extract mol related records and export for review, alternatively export the file for exclusion review
      if (nrow(object$mortality_filtered %>% select(conc1_unit) %>% filter(conc1_unit %like% "mol/L")) > 0){
        object <- export_mol_units(object = object)
        object$state <- 2
        project$object <- object
        save_project(project, project_path, save_project_steps)
        return(project)
      } else {
        #object <- export_exclude_list(object,project$project_path)
        export_chemical_list(object, project$project_path)
        object$state <- 3
        project$object <- object
        save_project(project, project_path, save_project_steps)
        return(project)
        }

    }

    # workflow step 3

    if(object$state == 2){
      # update records with mol related concentrations
      if (file.exists(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mol_weight.csv")))){
        object <- update_mol_units(object = object)
        export_chemical_list(object, project$project_path)
        object$state <- 3
        project$object <- object
        save_project(project, project_path, save_project_steps)
        return(project)
      }
    }

    if(object$state == 3){
      object <- remove_excluded_chemicals(object, project$project_path)
      if(update_chemicals == TRUE){
        object <- update_chemical_list(object, project_path = project_path, database_path = database_path)
      }
      object <- build_final_list(object)
      object <- calculate_water_solubility(object)
      object <- calculate_hours(object)

      #object <- calculate_water_solubility(object)
      object$state <- 4
      fwrite(object$results,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(object$parameters$ecotox_group), "_final_results.csv")))), sep = ",", dec = ".", na = NA)
      project$object <- object
      save_project(project, project_path,save_project_steps)
      message("[EcoToxR]:  The data pre-processing is finalised.")
      message("[EcoToxR]:  Please run the pivot export workflow.")
      return(project)
    }


}


