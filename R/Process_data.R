# functions

  process_data <- function(project,
                           all_species = FALSE,
                           dosing_group = "water_concentration",
                           duration_d = NA,
                           duration_h = NA,
                           duration_m = NA,
                           ecotoxgroup = NA,
                           effects = NA, #
                           habitat = NA, #
                           kingdoms = NA,
                           measurements = NA,
                           max_d = NA,
                           min_d = NA,
                           max_h = NA,
                           min_h = NA,
                           min_m = NA,
                           max_m = NA,
                           remove_operator_records = FALSE,
                           remove_formulation = FALSE,
                           sample_size = NA, # required for testing only
                           save_project_steps = TRUE,
                           species_selection = NA,
                           state = NA,
                           update_chemicals = FALSE,
                           water_solubility_cut_off = FALSE,
                           water_solubility_threshold = NA
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
        object$parameters$remove_operator_records <- remove_operator_records
        object$parameters$sample_size <- sample_size
        object$parameters$species_selection <- species_selection
        object$parameters$update_chemicals <- update_chemicals
        object$parameters$water_solubility_cut_off <- water_solubility_cut_off
        object$parameters$water_solubility_threshold <- water_solubility_threshold
    }

    # to be integrated
    #

    # Leftover from prepare data, to be considered or deleted later
    # object$all_selected_effects <-  object$all_selected_effects %>% drop_na(endpoint)
    #
    # object$all_selected_effects <- data.table(reduce(list(
    #     object$chemicals[, colnames(object$chemicals) %in% c("cas_number", "chemical_name", "dtxsid"), with = FALSE],
    #     object$all_selected_effects
    # ),
    # left_join,
    # by = "cas_number"))
    #
    # colnames(object$all_selected_effects)[2] <- "chemical_name"


    # Workflow step 1
    if (is.null(object$state) == TRUE) {

      # select only water organisms and merge tables
      message("[EcoToxR]:  Subsetting to habitat.")
      object$results_filtered <- object$merged_results %>% filter(organism_habitat %in% habitat)

      # Do preliminary filtering (group and effects)
      if (!is.na(ecotoxgroup)) {
          message("[EcoToxR]:  Filtering the datasets by the ecotox group (species group).")
          object$results_filtered <- object$results_filtered %>% filter(ecotox_group %like% ecotoxgroup)
      }

      # Too complicated, needs refactorization see issue #9 https://git.ufz.de/WANA/REcoTox/-/issues/9
      if (is.na(dosing_group)) {
          message("[EcoToxR]:  The dosing group cannot be empty or set to \"NA\"")
      } else if (dosing_group == "water_concentration"){
          message("[EcoToxR]:  Reducing to mass/water related units.")
          unit_pattern <- "(^(|(A|a)(e|i|E|I) )(|m|n|u|p)(g|mol)/(|m|u|p|d|)(L|l|m3)$)|(^(p)(p)(m|t|b)$)"
          object$results_filtered <- object$results_filtered %>% filter(conc1_unit %like% unit_pattern)

          message("[EcoToxR]:  The following units have been found in the data:")
          print(pull(object$results_filtered %>% select(conc1_unit) %>% unique() %>% group_by(conc1_unit) %>% arrange(.by_group = TRUE)))
      } else if (!dosing_group == "water_concentration") {
          message("[EcoToxR]:  The dosing group must be set to \"water_concentration\"")
      }


      # Subset algae data to kingdoms
      if (ecotoxgroup == "Algae" && !is.na(kingdoms)) {
          object$results_filtered <- object$results_filtered %>% filter(kingdom %in% kingdoms)
      }


      # Select by effect
      if (!is.logical(effects)) {
          message("[EcoToxR]:  Filtering the datasets by the effects")
          object$results_filtered <- object$results_filtered %>% filter(effect %in% effects)

          message("[EcoToxR]:  The following effects are available:")
          print(pull(object$results_filtered %>% select(effect) %>% unique() %>% group_by(effect) %>% arrange(.by_group = TRUE)))
      }

      # Select by measurement
      if (!is.logical(measurements)) {
          message("[EcoToxR]:  Filtering the datasets by the measurement.")
          object$results_filtered <- object$results_filtered %>% filter(measurement %in% measurements)

          message("[EcoToxR]:  The following measurements are available:")
          print(pull(object$results_filtered %>% select(measurement) %>% unique() %>% group_by(measurement) %>% arrange(.by_group = TRUE)))
      }

      # Remove formulations
      if (remove_formulation == TRUE) {
          message("[EcoToxR]:  Removing formulations.")
          object$results_filtered <- object$results_filtered %>% filter(conc1_type == "A")
          removed_formulation <- object$results_filtered %>% filter(conc1_type == "T")
          message(paste0("[EcoToxR]:  The dataset includes ", nrow(object$results_filtered), " records of active ingredients."))
          message(paste0("[EcoToxR]:  ", nrow(removed_formulation), " record(s) of formulations was/were removed."))
          rm(removed_formulation)
      }



    # Handling durations
      # Filter the time groups (minutes, hours, days)
    results_per_days <- object$results_filtered %>% filter(obs_duration_unit %in% duration_d)
    results_per_hours <- object$results_filtered %>% filter(obs_duration_unit %in% duration_h)
    results_per_minutes <- object$results_filtered %>% filter(obs_duration_unit %in% duration_m)


    # Subset to minimal and maximal durations

    results_per_days <- results_per_days %>% filter(obs_duration_mean >= min_d)
    results_per_days <- results_per_days %>% filter(obs_duration_mean <= max_d)

    results_per_hours <- results_per_hours %>% filter(obs_duration_mean >= min_h)
    results_per_hours <- results_per_hours %>% filter(obs_duration_mean <= max_h)

    results_per_minutes <- results_per_minutes %>% filter(obs_duration_mean >= min_m)
    results_per_minutes <- results_per_minutes %>% filter(obs_duration_mean <= max_m)

    object$results_filtered <- bind_rows(results_per_days, results_per_hours, results_per_minutes)

    rm(results_per_days, results_per_hours, results_per_minutes)


    # Prepare and export a list of species
    # species_selection = c("all", "selected", "standard_test_species")
    if (is.na(species_selection) | species_selection == "all") {
        species <- object$results_filtered %>% select(species_number, ecotox_group) %>%
            unique() %>%
            select(species_number) %>%
            left_join(object$species, by = "species_number")

        species <- species %>%
            add_column("include_species" = 1, .before = 1) %>%
            add_column("count_of_records" = 0, .after = 1)

        for (i in 1:nrow(species)) {
            species$count_of_records[i] <- length(which(object$results_filtered[, "species_number"] ==
                                                            species[i, "species_number"][[1]]))
        }
        species <- species %>% arrange(-count_of_records)

        message("[EcoToxR]:  The table with the species was written to the project folder.")
        #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
        write_csv(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))))
    } else if (species_selection == "manual") {

      species <- object$results_filtered %>%
          select(species_number, ecotox_group) %>%
          unique() %>%
          select(species_number) %>%
          left_join(object$species, by = "species_number") %>%
          add_column("include_species" = 0, .before = 1) %>%
          add_column("count_of_records" = 0, .after = 1)


          for (i in 1:nrow(species)) {
            species$count_of_records[i] <- length(which(object$results_filtered[, "species_number"] ==
                                                                    species[i, "species_number"][[1]]))
          }

          species <- species %>% arrange(-count_of_records)

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          write_csv(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      } else if (species_selection == "standard_test_species") {

          species <- object$results_filtered %>% select(species_number, ecotox_group) %>%
              filter(ecotox_group %like% "Standard Test Species") %>%
              unique() %>%
              select(species_number) %>%
              left_join(object$species, by = "species_number")

          species <- species %>%
              add_column("include_species" = 1, .before = 1) %>%
              add_column("count_of_records" = 0, .after = 1)

          for (i in 1:nrow(species)) {
              species$count_of_records[i] <- length(which(object$results_filtered[, "species_number"] ==
                                                              species[i, "species_number"][[1]]))
          }
          species <- species %>% arrange(-count_of_records)

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          write_csv(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      } else {
          species <- object$results_filtered %>% select(species_number, ecotox_group) %>%
              unique() %>%
              select(species_number) %>%
              left_join(object$species, by = "species_number")

          species <- species %>%
              add_column("include_species" = 1, .before = 1) %>%
              add_column("count_of_records" = 0, .after = 1)

          for (i in 1:nrow(species)) {
              species$count_of_records[i] <- length(which(object$results_filtered[, "species_number"] ==
                                                              species[i, "species_number"][[1]]))
          }
          species <- species %>% arrange(-count_of_records)

          message("[EcoToxR]:  Species selection has unvalid input. The table with the species was written to the project folder for manual review.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          write_csv(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotoxgroup), "_species_selection.csv")))))

      }

    # Select ecotox/species group (e.g. "Algae")
    message("[EcoToxR]:  The following ecotox/species group(s) are included:")
    print(pull(species %>% select(ecotox_group) %>% unique() %>% group_by() %>% arrange(.by_group = TRUE)))


      endpoints <- object$results_filtered %>%
          select(endpoint) %>%
          unique() %>%
          add_column("include_endpoint" = 1, .before = 1) %>%
          add_column("count_of_records" = 0, .after = 1)

      for (i in 1:nrow(endpoints)) {
          endpoints$count_of_records[i] <- length(which(object$results_filtered[, "endpoint"] ==
                                                            endpoints[i, "endpoint"][[1]]))
      }

      endpoints <- endpoints %>% group_by(endpoint) %>% arrange(.by_group = TRUE)

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
    endpoints <- read_csv(file.path(project_path, paste0(tolower(ecotoxgroup), "_endpoint_selection.csv")), show_col_types = FALSE)

    if(!all(unique(endpoints$include_endpoint) %in% c(0, 1))) {
        stop("[EcoToxR]:  Only 0 or 1 is allowed in the endpoint selection filter. Check the endpoint selection list and re-run the workflow.")
    }

    #if(unique(endpoints$include_endpoint)[[1]] == 0) {
    #    stop("[EcoToxR]:  There are only 0 in the endpoint selection filter. Check the species selection list and re-run the workflow.")
    #}

    object$results_filtered <- object$results_filtered %>%
      left_join(endpoints %>% select(include_endpoint, endpoint), by = "endpoint") %>%
      filter(include_endpoint == 1)

    # Filter the selected species

    species <- read_csv(file.path(project_path, paste0(tolower(ecotoxgroup), "_species_selection.csv")), show_col_types = FALSE)

    #if(!all(unique(species$include_species) %in% c(0, 1))) {
    #  stop("[EcoToxR]:  Only 0 or 1 is allowed in the species selection filter. Check the species selection list and re-run the workflow.")
    #}

    #if(unique(species$include_species) == 0) {
    #  stop("[EcoToxR]:  There are only 0 in the species selection filter. Check the species selection list and re-run the workflow.")
    #}

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
        object <- update_mol_units(object = object, database_path = project$database_path, project_path = project$project_path)
        export_chemical_list(object, project$project_path)
        object$state <- 3
        project$object <- object
        save_project(project, save_project_steps)
        return(project)
      }
      #else {
          # needs an alternative procedure, e.g. check for mol/L and update stage
      #}
    }

    if(object$state == 3){

      # This routine just exports an exclusion and inclusion list, but does not remove chemicals from chemical list
      object <- remove_excluded_chemicals(object, project$project_path)

      if(update_chemicals == TRUE){
        object <- update_chemical_list(object, project_path = project_path, database_path = database_path)
      }
      object <- build_final_list(object)
      object <- calculate_hours(object)

      #object <- calculate_water_solubility(object)
      object$state <- 4
      write_csv(x = object$results, file = suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_final_results.csv")))))
      project$object <- object
      save_project(project, save_project_steps)
      message("[EcoToxR]:  The data pre-processing is finalised.")
      message("[EcoToxR]:  Please run the pivot export workflow.")
      return(project)
    }
  }

