# helper functions

build_final_list <- function(object = object){

  results <- data.table(object$mortality_filtered)

  if(is.null(results$include_species)){
    results[, "include_species" := NA]
  }

  select_columns <- c("cas_number", "cas", "chemical_name", "compound_class", "test_location", "conc1_type", "conc1_mean_op",
                      "conc1_mean", "conc1_min_op", "conc1_min", "conc1_max_op", "conc1_max", "conc1_unit",
                      "conc1_comments", "concentration_mean", "concentration_unit", "obs_duration_mean",
                      "obs_duration_unit", "test_id", "reference_number", "endpoint", "endpoint_comments", "effect",
                      "effect_comments", "measurement", "measurement_comments", "organism_lifestage", "common_name", "latin_name", "kingdom", "phylum_division",
                      "subphylum_div", "superclass", "class", "tax_order", "family", "genus", "species", "subspecies",
                      "variety", "ecotox_group", "result_id", "reference_db", "reference_type", "author",
                      "title", "source", "publication_year", "include_endpoint", "include_species")

  results <- data.table(results[, ..select_columns])

  results <- data.table(left_join(results, object$chemprop[, c("cas_number", "DTXSID", "PREFERRED_NAME", "CASRN", "SMILES",
                                                  "INCHIKEY", "AVERAGE_MASS", "OPERA_LOG_P", "OPERA_LOG_P_AD",
                                                  "OPERA_LOG_D_74", "OPERA_LOG_D_AD",	"OPERA_LOG_S_74",	"OPERA_LOG_S_AD",
                                                  "ACD_LOG_P", "ACD_LOG_D_74", "ACD_LOG_S_74", "JC_LOG_P",	"JC_LOG_D_74",
                                                  "JC_LOG_S_74", "EXCLUDE", "REMARKS")],
                                                  by = "cas_number")
  )

  # reorder list

  select_columns <- c("EXCLUDE", "REMARKS", "cas_number", "cas", "DTXSID", "PREFERRED_NAME", "test_location", "reference_number",
                      "conc1_type", "conc1_mean_op", "conc1_mean", "conc1_min_op", "conc1_min", "conc1_max_op", "conc1_max",
                      "conc1_unit", "conc1_comments", "concentration_mean", "concentration_unit", "test_id",
                      "endpoint", "endpoint_comments", "effect", "effect_comments",
                      "measurement", "measurement_comments", "obs_duration_mean",
                      "obs_duration_unit", "organism_lifestage", "common_name", "latin_name", "kingdom",
                      "phylum_division", "subphylum_div", "superclass", "class", "tax_order", "family",
                      "genus", "species", "subspecies", "variety", "ecotox_group", "chemical_name",
                      "compound_class", "CASRN", "SMILES", "INCHIKEY", "AVERAGE_MASS",
                      "OPERA_LOG_P", "OPERA_LOG_P_AD", "OPERA_LOG_D_74", "OPERA_LOG_D_AD", "OPERA_LOG_S_74", "OPERA_LOG_S_AD",
                      "ACD_LOG_P", "ACD_LOG_D_74", "ACD_LOG_S_74", "JC_LOG_P", "JC_LOG_D_74",
                      "JC_LOG_S_74", "reference_db", "result_id", "reference_type", "author",
                      "title", "source", "publication_year", "include_endpoint", "include_species")


  results <- results[,..select_columns]

  object$results <- data.table(results)
  return(object)
}


calculate_hours <- function(object = object){
  message("[EcoToxR]:  Recalculating duration data (day to hours).")
  results <- data.table(object$results)

  results_h <- results %>% filter(obs_duration_unit == "h")
  results_d <- results %>% filter(obs_duration_unit == "d")

  results_d <- results_d %>% filter(obs_duration_unit == "d") %>% mutate(obs_duration_mean = 24 * obs_duration_mean) %>% mutate(obs_duration_unit = "h")

  results <- data.table(rbind(results_h, results_d))

  object$results <- results

  return(object)

}


calculate_water_solubility <- function(object = object){
  results <- data.table(object$results)

  # OPERA
  results <- results %>% mutate(OPERA_S_mg_L = 10^OPERA_LOG_S_74 * AVERAGE_MASS * 1000) %>% filter(is.na(EXCLUDE)) # OPERA output g/L

  # ACD
  results <- results %>% mutate(ACD_S_mg_L = 10^ACD_LOG_S_74 * AVERAGE_MASS * 1000) %>% filter(is.na(EXCLUDE)) # ACD output g/L

  # JC
  results <- results %>% mutate(JC_S_mg_L = 10^JC_LOG_S_74 * AVERAGE_MASS*1000) %>% filter(is.na(EXCLUDE)) # JC output g/L

  # Domain estimate solubility domain (based on ideas in ChemProp)
  # Case 1: if EC <= Sw -> 3
  # Case 2: if 5*log10 Sw >= EC > Sw -> 2
  # Case 3: if 10*log10 Sw >= EC > 5*log10 Sw -> 1
  # Case 4: if EC > 10*log10 Sw -> 0

  s <- c("OPERA_S_mg_L", "ACD_S_mg_L", "JC_S_mg_L")

  length_progressbar <- length(s)
  pb <- progress::progress_bar$new(
    format = "[EcoToxR]:  Estimating the solubility domain [:bar] :percent ETA: :eta",
    total = length_progressbar, clear = FALSE, width = 80)

  for(i in s){

  AD <- paste0(i, "_AD")
  results[, "AD"] <- NA


  # If concentration_mean is <= Sw, the
  for(j in 1:nrow(results)){
    if(is.na(results$EXCLUDE[j]) & is.na(results[j, ..i])){
      results[j, AD] <- NA
      }

    else if(is.na(results$EXCLUDE[j]) & results$concentration_mean[j] <= results[j, ..i]){
      results[j, AD] <- 3
      }

    else if(results$concentration_mean[j] > results[j, ..i] & results$concentration_mean[j] <= 10^5 * log10(results[j, ..i])){
      results[j, AD] <- 2
      }

    else if (results$concentration_mean[j] > 10^5 * log10(results[j, ..i]) & results$concentration_mean[j] <= 10^10 * log10(results[j, ..i])){
      results[j, AD] <- 1
      }

    else if (results$concentration_mean[j] > 10^10 * log10(results[j, ..i])){
      results[j, AD] <- 0
    }
    }
  pb$tick()
  }

  object$results <- data.table(results)
  return(object)
}


export_chemical_list <- function(object,project_path){
  message("[EcoToxR]:  Exporting the list of included chemicals to the project folder.")
  message("[EcoToxR]:  Check the file for exclusion of chemicals.")
  message("[EcoToxR]:  Impute missing values for physical-chemical properties.")
  message(paste0("[EcoToxR]:  Edit the file ",tolower(object$parameters$ecotoxgroup),"_chemical_list.csv and re-run the workflow."))
  chemical_list <- data.table(unique(object$mortality_filtered[,c("cas_number", "cas", "chemical_name")]))
  chemical_list <- left_join(chemical_list,object$chemprop[,c("cas_number", "FOUND_BY", "DTXSID", "PREFERRED_NAME",
                                                               "CASRN", "INCHIKEY", "IUPAC_NAME", "SMILES", "INCHI_STRING",
                                                               "MOLECULAR_FORMULA", "AVERAGE_MASS", "MONOISOTOPIC_MASS",
                                                               "MS_READY_SMILES", "QSAR_READY_SMILES", "OPERA_LOG_P", "OPERA_LOG_P_AD",
                                                              "OPERA_LOG_D_74", "OPERA_LOG_D_AD", "OPERA_LOG_S_74", "OPERA_LOG_S_AD",
                                                              "ACD_LOG_P", "ACD_LOG_D_74", "ACD_LOG_S_74",
                                                              "JC_LOG_P", "JC_LOG_S_74", "JC_LOG_D_74", "EXCLUDE", "REMARKS")],
                                                               by = "cas_number"
                             )

  chemical_list <- chemical_list[order(chemical_list$chemical_name), ]
  fwrite(chemical_list, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")))), sep = ",", dec = ".", na = NA)

}

export_exclude_list <- function(object,project_path){
  message("[EcoToxR]:  Exporting the final list to the project folder.")
  message("[EcoToxR]:  Check the file for exclusion of records.")
  message(paste0("[EcoToxR]:  Edit the file ",tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv and re-run the workflow."))
  exclude <- data.table(matrix(nrow = nrow(object$mortality_filtered), ncol = 2))
  colnames(exclude) <- c("exclude", "exclusion_comment")
  exclude_list <- data.table(exclude,object$mortality_filtered)
  object$mortality_filtered_exclude_list <- data.table(exclude_list[, c("exclude", "exclusion_comment", "result_id", "cas_number", "cas", "chemical_name",
                                                                       colnames(exclude_list)[grep(pattern = "conc1",colnames(exclude_list))],
                                                                       "species_number", "concentration_mean", "concentration_unit",
                                                                       "latin_name", "author", "title", "source", "publication_year"),
                                                                    with = FALSE])
  object$mortality_filtered_exclude_list <- data.table(object$mortality_filtered_exclude_list[order(chemical_name)])
  fwrite(data.table(unique(object$mortality_filtered_exclude_list$cas)), suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list_cas.csv")))), sep = ",", dec = ".")
  fwrite(object$mortality_filtered_exclude_list, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv")))), sep = ",", dec = ".")
  return(object)
}

remove_excluded_chemicals <- function(object, project_path){

  chemical_list <- data.table(fread(file.path(project_path,paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")),
                                    sep = ",", dec = ".", na.strings = c("NA", "", NA, NaN)))

  suppressWarnings(chemical_list <- format_chemical_properties(chemical_list))
  object$chemical_list <- chemical_list


  if(!all(unique(chemical_list$exclude) == 1, na.rm = TRUE)){
    stop("Only NA or 1 is allowed in the exclusion filter. Check the exclusion list and re-run the workflow.")
  }
  inclusion_list <- chemical_list[, chemical_list[is.na(EXCLUDE)]]$cas_number
  exclusion_list <- chemical_list[chemical_list[, EXCLUDE %like% 1]]$cas_number

  object$mortality_removed_chemicals <- object$mortality_filtered[cas_number %in% exclusion_list]
  object$mortality_filtered <- object$mortality_filtered[cas_number %in% inclusion_list]

  fwrite(object$mortality_filtered,file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_processed.csv")),
                          sep = ",", dec = ".")

  if(length(exclusion_list) > 0){
    fwrite(object$mortality_filtered, suppressWarnings(normalizePath(file.path(project_path,paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_processed_excluded.csv")))), sep = ",", dec = ".")
  }
  return(object)
}

update_chemical_list <- function(object = object, project_path = project_path, database_path = database_path){

  # some kind dirty solution...
  #object$chemprop <- object$chemprop[chemprop, on = c("cas_number"), AVERAGE_MASS := i.AVERAGE_MASS]
  #object$chemprop <- object$chemprop %>% mutate_all(na_if,"")

  chemprop <- data.table(object$chemprop)
  chemical_list <- data.table(object$chemical_list)
  chemical_list <- format_chemical_properties(chemical_list)
  chemprop <- data.table(merge(chemprop, chemical_list, by='cas_number', all = TRUE, suffixes = c("", ".update")))
  col_names1 <- grep(".update", colnames(chemprop))
  col_names2 <- colnames(chemprop)[grep(".update", colnames(chemprop))]
  col_names3 <- gsub(".update", replacement = "", col_names2)



  length_progressbar <- length(col_names1)
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Updating the chemical properties [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)


  for(i in 1:length(col_names1)){
    c1 <- which(colnames(chemprop) == col_names3[i])
    c2 <- which(colnames(chemprop) == col_names2[i])
    for(j in 1:nrow(chemprop)){
      if(!is.na(chemprop[j,..c2])){
      chemprop[j, c1] <- chemprop[j, ..c2]
      } else next()
    }
    pb$tick()
  }

chemprop <- data.table(chemprop)
chemprop[,(col_names2) := NULL]
export_chemical_properties(object = chemprop, database_path = database_path, project_path = project_path)
object$chemprop <- chemprop

return(object)
}

remove_asterics <- function(object){
  # message("Remove effect reportings which have asterics")
  # Remove all wired effect reportings
  ## Replace by an export mechanism
  nrow(object)

  suppressWarnings(
    {
      object$conc1_mean <- as.numeric(object[, conc1_mean])
      object$conc1_min <- as.numeric(object[, conc1_min])
      object$conc1_max <- as.numeric(object[, conc1_max])
    }
  )

  conc_remove <- c("[*]")
  object <- object[!conc1_mean %like% conc_remove]
  object <- object[!conc1_min %like% conc_remove]
  object <- object[!conc1_max %like% conc_remove]
  #object <- object[!conc1_min_op %like% conc_remove]
  #object <- object[!conc1_max_op %like% conc_remove]
  #object$conc_remove <- c("~", ">=")
  #object$all_selected_effects <- object$all_selected_effects[!conc1_mean_op %in% object$conc_remove]

  object <- filter(object, !is.na(conc1_mean))

  return(object)
}


# Handle the chemical lists
export_chemical_properties <- function(object, database_path = database_path,
                                       project_path = project_path){

  object <- format_chemical_properties(object)
  fwrite(object,suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))), sep = ",", dec = ".", na = NA)
  fwrite(object,suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))), sep = ",", dec = ".", na = NA)

}

create_chemical_properties <- function(database_path){
  if (file.exists(file.path(database_path,"chemical_properties.csv"))){
    message("[EcoToxR]:  Reading chemical properties (this is a custom file).")
    object <- fread(file.path(database_path,"chemical_properties.csv"), sep = ",", dec = ".", na.strings = c("NA", "", NA, NaN))
    object <- format_chemical_properties(object)

  } else {

    message("[EcoToxR]:  A custom file for the storage of chemical properties (mol weight) is compiled.")
    message("[EcoToxR]:  The file is stored in the current database folder and will be copied to the project folder for backup.")
    object <- data.table(matrix(nrow = 0, ncol = 32))
    colnames(object) <- c("cas_number", "cas", "chemical_name", "FOUND_BY", "DTXSID", "PREFERRED_NAME",
                          "CASRN", "INCHIKEY", "IUPAC_NAME", "SMILES", "INCHI_STRING", "MOLECULAR_FORMULA",
                          "AVERAGE_MASS", "MONOISOTOPIC_MASS", "MS_READY_SMILES", "QSAR_READY_SMILES",
                          "OPERA_LOG_P", "OPERA_LOG_P_AD", "OPERA_LOG_D_74",
                          "OPERA_LOG_D_AD", "OPERA_LOG_S_74", "OPERA_LOG_S_AD",
                          "ACD_LOG_P", "ACD_LOG_D_74", "ACD_LOG_S_74",
                          "JC_LOG_P", "JC_LOG_S_74", "JC_LOG_D_74", "EXCLUDE", "REMARKS"
                          )
    object <- format_chemical_properties(object)
    suppressWarnings(
      fwrite(object,suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))), sep = ",", dec = ".", na = NA)
    )

  }
  return(object)
}

format_chemical_properties <- function(object){
  object$cas_number <- as.integer(object$cas_number)
  object$cas <- as.character(object$cas)
  object$chemical_name <- as.character(object$chemical_name)
  object$FOUND_BY <-  as.character(object$FOUND_BY)
  object$DTXSID <- as.character(object$DTXSID)
  object$PREFERRED_NAME <- as.character(object$PREFERRED_NAME)
  object$CASRN <- as.character(object$CASRN)
  object$INCHIKEY <- as.character(object$INCHIKEY)
  object$IUPAC_NAME <- as.character(object$IUPAC_NAME)
  object$SMILES <- as.character(object$SMILES)
  object$INCHI_STRING <- as.character(object$INCHI_STRING)
  object$MOLECULAR_FORMULA <- as.character(object$MOLECULAR_FORMULA)
  object$AVERAGE_MASS <- as.numeric(object$AVERAGE_MASS)
  object$MONOISOTOPIC_MASS <- as.numeric(object$MONOISOTOPIC_MASS)
  object$MS_READY_SMILES <- as.character(object$MS_READY_SMILES)
  object$QSAR_READY_SMILES <- as.character(object$QSAR_READY_SMILES)
  object$OPERA_LOG_P <- as.numeric(object$OPERA_LOG_P)
  object$OPERA_LOG_P_AD <- as.numeric(object$OPERA_LOG_P_AD)
  object$OPERA_LOG_D_74 <- as.numeric(object$OPERA_LOG_D_74)
  object$OPERA_LOG_D_AD <- as.numeric(object$OPERA_LOG_D_AD)
  object$OPERA_LOG_S_74 <- as.numeric(object$OPERA_LOG_S_74)
  object$OPERA_LOG_S_AD <- as.numeric(object$OPERA_LOG_S_AD)
  object$ACD_LOG_P <- as.numeric(object$ACD_LOG_P)
  object$ACD_LOG_D_74 <- as.numeric(object$ACD_LOG_D_74)
  object$ACD_LOG_S_74 <- as.numeric(object$ACD_LOG_S_74)
  object$JC_LOG_P <- as.numeric(object$JC_LOG_P)
  object$JC_LOG_D_74 <- as.numeric(object$JC_LOG_D_74)
  object$JC_LOG_S_74 <- as.numeric(object$JC_LOG_S_74)
  object$EXCLUDE <- as.numeric(object$EXCLUDE)
  object$REMARKS <- as.character(object$REMARKS)
  return(object)
}


export_mol_units <- function(object, project_path = project$project_path){

  ecotoxgroup <- object$parameters$ecotoxgroup
  file_name <- suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))))

  message("[EcoToxR]:  Exporting the list of chemicals for mol/L to mg/L conversion to the project folder.")
  message(paste0("[EcoToxR]:  Edit the file ",tolower(ecotoxgroup),"_mol_weight.csv and re-run the workflow."))
  message("[EcoToxR]:  Add only missing data to column 'AVERAGE_MASS'.")

  chemprop <- unique(object$mortality_filtered[, c("cas_number", "cas", "chemical_name"), with = FALSE][which(object$mortality_filtered[, concentration_unit %like% "mol/L"])])
  chemprop$cas_number <- as.numeric(chemprop$cas_number)
  chemprop$cas <- as.character(chemprop$cas)
  chemprop$chemical_name <- as.character(chemprop$chemical_name)
  chemprop$AVERAGE_MASS <- NA

  if(nrow(object$chemprop)>0){
    suppressWarnings(chemprop <- left_join(chemprop, object$chemprop[, c(1, 4:21)], by = "cas_number"))
  }

  chemprop$AVERAGE_MASS.x <- NULL

  col_chemprop <- colnames(chemprop)
  col_num <- which(col_chemprop == "AVERAGE_MASS.y")
  colnames(chemprop)[col_num] <- "AVERAGE_MASS"


  chemprop <- chemprop[order(chemprop$chemical_name),]
  fwrite(chemprop,file_name, sep = ",", dec = ".", na = NA)
  return(object)
}

update_mol_units <- function(object, database_path = project$database_path, project_path = project$project_path){
  message("[EcoToxR]:  Recalulating record unit mol/L to mg/L.")
  ecotoxgroup <- object$parameters$ecotoxgroup
  file_name <- file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))



  chemprop <- data.table(fread(file_name, sep = ",", dec = ".", na.strings = c("NA", "", NA, NaN)))
  chemprop <- format_chemical_properties(chemprop)

  # subset
  nonmol_records <- object$mortality_filtered[!concentration_unit %like% "mol/L"]
  mol_records <- object$mortality_filtered[concentration_unit %like% "mol/L"]

  for(i in 1:nrow(mol_records)){
    cas_id <- chemprop[, which(cas_number %like% mol_records[i, cas_number])]

    mol_records[i, "concentration_mean"] <- mol_records[i, concentration_mean] * 1000 * chemprop[cas_id, AVERAGE_MASS]
    }

  object$mortality_filtered <- data.table(rbind(nonmol_records, mol_records))[order(cas_number)]


  object$mortality_filtered$concentration_unit[which(object$mortality_filtered[, concentration_unit %like% "mol/L"])] <- "mg/L"

  message("[EcoToxR]: The following units are ")
  print(unique(object$mortality_filtered$concentration_unit))

  #message("[EcoToxR]: Update the chemical properties.")

  object$chemprop <- object$chemprop[chemprop, on = c("cas_number"), AVERAGE_MASS := i.AVERAGE_MASS]

  object$chemprop <- object$chemprop %>% mutate_all(na_if, "")

  #fwrite(object$chemprop,file.path(database_path,"chemical_properties.csv"), sep = ",", dec = ".", quote = "\"", na = NA)
  fwrite(object$chemprop,suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))), sep = ",", dec = ".", na = c(NA, NaN, ""))

  return(object)
}

convert_units <- function(object, sample_size = NA){

  object <- data.table(object)

  if(!is.na(sample_size)){
    object <- object[sample(.N, sample_size)]
  }

  #message("Add new columns for concentrations")

  object$concentration_mean <- NA
  object$concentration_unit <- NA

  object$concentration_mean <- object$conc1_mean
  object$concentration_unit <- object$conc1_unit

  # update concentration units if NA


  length_progressbar <- length(object$concentration_mean[is.na(object$concentration_mean)])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Update missing concentration mean with mean from min and max [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:length(object$concentration_mean[is.na(object$concentration_mean)])){
      object[which(object[, is.na(concentration_unit)])][i]$concentration_mean <-
        mean(object[which(object[, is.na(concentration_unit)])][i]$conc1_min,object[which(object[, is.na(concentration_unit)])][i]$conc1_max)
      pb$tick()

    }
  }

  # convert similar unit to SI conform units (e.g. ppm to mg/L)

  # ng related
  message("[EcoToxR]:  Converting units to human readible format (i.e. mg/L).")
  message("[EcoToxR]:  The following units will be converted:")
  print(unique(object$concentration_unit))

  # pg/L
  pg <- "^pg/(l|L)$"
  length_progressbar <- nrow(object[which(object[, concentration_unit %like% pg])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting pg/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% pg])])){
      object[which(object[, concentration_unit %like% pg])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% pg])][i]$concentration_mean / 1e+09
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% pg])] <- "mg/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping pg/L like unit conversion")}


  #ng
  ng <- "^(|A(I|E|i|e) )ng/(l|L|dm3)$|^ppt$|^(|A(I|E|i|e) )pg/(ml|mL)$"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% ng])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting ng/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% ng])])){
      object[which(object[, concentration_unit %like% ng])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% ng])][i]$concentration_mean / 1e+06
      pb$tick()
    }


    object$concentration_unit[which(object[, concentration_unit %like% ng])] <- "mg/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping ng/L like unit conversion")}


  #µg related
  ug <- "(^(|(A|a)(I|E|i|e) )ug/(l|L|dm3)$)|^ppb$|(^pg/u(l|L)$)|(^(|(A|a)(I|E|i|e) )ng/m(l|L)$)|^pg/u(l|L)"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% ug])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting ug/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% ug])])){
      object[which(object[, concentration_unit %like% ug])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% ug])][i]$concentration_mean / 1e+03
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% ug])] <- "mg/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping ug/L like units conversion")}

  # g/L
  g <- "(^(|(A|a)(e|i|E|I) )(g/(l|L)$))|(^mg/(ml|mL)$|(^g/dm3$))"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% g])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting g/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% g])])){
      object[which(object[, concentration_unit %like% g])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% g])][i]$concentration_mean*1e+03
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% g])] <- "mg/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping g/L like units conversion")}



  message("[EcoToxR]:  Converting mg/L like units")
  mg <- "(^(|((A|a)(I|E|i|e)) )mg/(l|L|dm3)$)|^ppm$|(^(|((A|a)(I|E|i|e)) )ug/m(l|L|m3)$)|^g/m3$"
  object$concentration_unit[which(object[, concentration_unit %like% mg])] <- "mg/L"




  #pmol related
  pmol <- "(^(|(A|a)(I|E|i|e) )pmol/(l|L)$)"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% pmol])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting pmol/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% pmol])])){
      object[which(object[, concentration_unit %like% pmol])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% pmol])][i]$concentration_mean / 1e+12
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% pmol])] <- "mol/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping pmol/L like unit conversion")}



  #nmol related
  nmol <- "(^(|(A|a)(I|E|i|e) )nmol/(l|L)$)"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% nmol])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting nmol/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% nmol])])){
      object[which(object[, concentration_unit %like% nmol])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% nmol])][i]$concentration_mean / 1e+09
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% nmol])] <- "mol/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping nmol/L coversion")}

  #umol related
  umol <- "(^(|(A|a)(I|E|i|e) )umol/(dm3|l|L)$)|^nmol/ml$|^mmol/m3$"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% umol])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Coverting umol/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[,concentration_unit %like% umol])])){
      object[which(object[, concentration_unit %like% umol])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% umol])][i]$concentration_mean / 1e+06
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% umol])] <- "mol/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping umol/L like unit conversion")}

  #mmol related
  mmol <- "^(|(A|a)(I|E|i|e) )mmol/(dm3|l|L)$"

  length_progressbar <- nrow(object[which(object[, concentration_unit %like% mmol])])
  if(length_progressbar > 0){
    pb <- progress::progress_bar$new(
      format = "[EcoToxR]:  Converting mmol/L like units [:bar] :percent ETA: :eta",
      total = length_progressbar, clear = FALSE, width = 80)

    for(i in 1:nrow(object[which(object[, concentration_unit %like% mmol])])){
      object[which(object[, concentration_unit %like% mmol])][i]$concentration_mean <-
        object[which(object[, concentration_unit %like% mmol])][i]$concentration_mean / 1e+03
      pb$tick()

    }

    object$concentration_unit[which(object[, concentration_unit %like% mmol])] <- "mol/L"
    pb$terminate()
  } else {message("[EcoToxR]:  Skipping mmol/L like unit conversion")}

  message("[EcoToxR]:  The following units are left in the data:")
  print(unique(object$concentration_unit))


  return(object)

}

# clean_columns <- function(object){
#   na <- sapply(object, function(x) sum(length(which(is.na(x)))))
#   na <- names(na[na < nrow(object)])
#   object <- object[,na, with = FALSE]
#   return(object)
# }

# clean_rows <- function(object){
#   object[!is.na(object[,conc1_mean]),]
#   return(object)
# }

save_project <- function(object = object, project_path = project_path, save_project_steps = save_project_steps){
  state <- object$object$state
  ecotoxgroup <- object$object$parameters$ecotoxgroup
  if(isTRUE(save_project_steps)){
    project <- object
    message("[EcoToxR]:  Saving the project state.")
    save(project,file = file.path(project_path,paste0(tolower(ecotoxgroup),"_state", state,".RData")), compress = TRUE)
  }
}
