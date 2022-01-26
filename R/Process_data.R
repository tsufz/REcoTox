# functions

  process_data <- function(project,
                           ecotox_group = NA,
                           sample_size = NA,
                           all_species = FALSE,
                           measurements = c("MORT","GGRO","SURV"),
                           min_h = 0, min_d = 0,
                           max_h = 120, max_d = 5,
                           duration_h = c("h","ht","hph","hpf","hv","hbf"),
                           duration_d = c("d","dph","dpf"),
                           kingdoms = c("Chromista","Plantae","Monera"),
                           remove_asterics_values = TRUE,
                           water_solubility_cut_off = FALSE,
                           water_solubility_threshold = NA,
                           remove_operator_records = FALSE,
                           save_project_steps = TRUE,
                           update_chemicals = FALSE
                           ){
    object <- project$object


    # Record the settings
    object$parameters$ecotox_group <- ecotox_group
    object$parameters$measurements <- measurements
    object$parameters$min_h <- min_d
    object$parameters$min_d <- min_d
    object$parameters$max_h <- max_d
    object$parameters$max_d <- max_h
    object$parameters$duration_h <- duration_h
    object$parameters$duration_d <- duration_d
    object$parameters$kingdom <- kingdoms
    object$parameters$sample_size <- sample_size
    object$parameters$remove_asterics_values <- remove_asterics_values
    object$parameters$all_species <- all_species

    # Workflow step 1
    if(is.null(object$state) == TRUE){
      # Do preliminary filtering (group and effects)
      if (!is.na(object$parameters$ecotox_group <- ecotox_group)) {
          message("[EcoToxR]:  Filtering the datasets by the ecotox group (species group).")
          all_selected_effects <- data.table(object$all_selected_effects[ecotox_group %like% ecotox_group])
       }

      # Select mortality records only
      message("[EcoToxR]:  The following ecotox group(s) are included:")
      print(unique(all_selected_effects$ecotox_group))

      message("[EcoToxR]:  Filtering the datasets by the measurement.")
      object$mortality <- data.table(all_selected_effects %>% filter(measurement %in% measurements))
      message("[EcoToxR]:  The following measurements are included:")
      print(sort(unique(object$mortality$measurement)))


      # Subset algae data to kingdom

      if(ecotox_group == "Algae"){
        object$mortality <- object$mortality[kingdom %in% kingdoms]
      }


      # Select given durations
      mortality_h <- object$mortality[obs_duration_unit %in% duration_h]
      mortality_d <- object$mortality[obs_duration_unit %in% duration_d]

      # Subset to minimal and maximal durations

      mortality_h <- subset(mortality_h, obs_duration_mean >= min_h)
      mortality_h <- subset(mortality_h, obs_duration_mean <= max_h)

      mortality_d <- subset(mortality_d, obs_duration_mean >= min_d)
      mortality_d <- subset(mortality_d, obs_duration_mean <= max_d)

      object$mortality <- data.table(rbind(mortality_h, mortality_d))

      mortality_h <- NULL
      mortality_d <- NULL

      # Prepare and export a list of endpoints
      endpoints <- data.table(unique(object$mortality$endpoint))

      # Prepare and export a list of species
      if(!all_species){

          species <- data.table(unique(object$mortality$species_number))
          colnames(species)[1] <- "species_number"

          species <- data.table(reduce(list(
            species,
            object$species
          ),
          left_join,
          by = "species_number"))

          include <- data.table(matrix(nrow = nrow(species), ncol = 2))
          colnames(include) <- c("include_species","count_of_records")
          species <- data.table(include, species)
          species <- species[order(latin_name),]
          include <- NULL

          for(i in 1:nrow(species)){
            species$count_of_records[i] <- length(which(object$mortality[, species_number] ==
                                                                    species[, species_number][i]))
          }

          message("[EcoToxR]:  The table with the species was written to the project folder.")
          #message("[EcoToxR]:  Please edit the file and rerun the workflow.")
          fwrite(species, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group),"_species_selection.csv")))),  sep = ",", dec = ".")
      }


      endpoints <- data.table(unique(object$mortality$endpoint))
      include <- data.table(matrix(nrow = nrow(endpoints), ncol = 1))
      colnames(include) <- "include_endpoint"
      colnames(endpoints)[1] <- "endpoint"
      endpoints <- data.table(include, endpoints)
      endpoints <- endpoints[order(endpoint), ]
      include <- NULL


      for(i in 1:nrow(endpoints)){
        endpoints$count[i] <- length(which(object$mortality[, endpoint] ==  endpoints[, endpoint][i]))

      }

      message("[EcoToxR]:  The table with the endpoints was written to the project folder.")
      message("[EcoToxR]:  Edit the file(s) and re-run the workflow.")
      fwrite(endpoints, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group), "_endpoint_selection.csv")))),  sep = ",", dec = ".")

      object$state <- 1

      project$object <- object
      save_project(project, project_path, save_project_steps)

      fwrite(object$mortality,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group), "_mortality_unfiltered.csv")))),  sep = ",", dec = ".")
      return(project)
    }

    # Workflow step 2

    if(object$state == 1){
      endpoints <- data.table(fread(file.path(project_path,paste0(tolower(ecotox_group), "_endpoint_selection.csv"))))
      if(!all(unique(endpoints$include) == 1, na.rm = TRUE)){
        stop("[EcoToxR]:  Only NA or 1 is allowed in the endpoint selection filter. Check the endpoint selection list and re-run the workflow.")
      }

      if(isTRUE(is.na(unique(endpoints$include)))){
        stop("[EcoToxR]:  There are only NA in the endpoint selection filter. Check the species selection list and re-run the workflow.")
      }

      object$mortality <- left_join(object$mortality, endpoints[, c(1, 2)], by = "endpoint") # merge
      object$mortality_filtered <- data.table((filter(object$mortality, include_endpoint == 1)))
      object$mortality_filtered_removed_endpoint <- filter(object$mortality, is.na(include_endpoint))

      # Remove asterics in the concentration fields
      if(remove_asterics_values == TRUE){
        object$mortality_filtered <- remove_asterics(object$mortality_filtered)
      }

      if(all_species == FALSE){
        species <- data.table(fread(file.path(project_path,paste0(tolower(ecotox_group), "_species_selection.csv")), sep = ",", dec = "."))
        species <- species[order(include_species),]
        if(!all(unique(species$include) %in% c(1, NA))){
          stop("[EcoToxR]:  Only NA or 1 is allowed in the species selection filter. Check the species selection list and re-run the workflow.")
        }
        if(is.na(unique(species$include)[1])){
          stop("[EcoToxR]:  There are only MA in the species selection filter. Check the species selection list and re-run the workflow.")
        }

        object$mortality_filtered <- left_join(object$mortality_filtered, species[, c(1, 3)], by = "species_number")
        object$mortality_filtered <- data.table(filter(object$mortality_filtered, include_species == 1))
        object$mortality_filtered_removed_species <- filter(object$mortality_filtered, is.na(include_species))
        fwrite(object$mortality_filtered,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group),"_mortality_filtered.csv")))),  sep = ",", dec = ".")
        fwrite(object$mortality_filtered,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group),"_mortality_removed_species.csv")))),  sep = ",", dec = ".")
      }

      # Check measurements

      # Convert the units and save the file
      object$mortality_filtered <- data.table(convert_units(object$mortality_filtered, sample_size))
      fwrite(object$mortality_filtered_removed_endpoint,suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(ecotox_group), "_mortality_removed_endpoints.csv")))),  sep = ",", dec = ".")

      # extract mol related records and export for review, alternatively export the file for exclusion review
      if (nrow(object$mortality_filtered[concentration_unit %like% "mol/L"]) > 0){
        object <- export_mol_units(object = object)
        object$state <- 2
        project$object <- object
        save_project(project, project_path, save_project_steps)
        return(project)
      } else {
        #object <- export_exclude_list(object,project$project_path)
        export_chemical_list(object,project$project_path)
        object$state <- 3
        project$object <- object
        save_project(project, project_path,save_project_steps)
        return(project)
        }

    }

    # workflow step 3

    if(object$state == 2){
      # update records with mol related concentrations
      if (file.exists(file.path(project_path, paste0(tolower(object$parameters$ecotox_group), "_mol_weight.csv")))){
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


