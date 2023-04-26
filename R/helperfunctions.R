# helper functions

build_final_list <- function(object = object){

  results <- object$results_filtered

  if(is.null(results$include_species)){
    results[, "include_species" := NA]
  }


  results <- results %>% select(cas_number, cas, chemical_name, compound_class,
                                dtxsid_ecotox, test_location, conc1_type, conc1_mean_op,
                                conc1_mean, conc1_min_op, conc1_min, conc1_max_op,
                                conc1_max, conc1_unit, conc1_comments,
                                concentration_mean, concentration_unit, obs_duration_mean,
                                obs_duration_unit, test_id, reference_number,
                                endpoint, endpoint_comments, effect,
                                effect_comments, measurement, measurement_comments,
                                organism_lifestage, common_name, latin_name, kingdom,
                                phylum_division, subphylum_div, superclass, class,
                                tax_order, family, genus, species, subspecies,
                                variety, ecotox_group, result_id, reference_db,
                                reference_type, author, title, source,
                                publication_year, include_endpoint, include_species) %>%
    rename(dtxsid_ecotox = dtxsid_ecotox)


  results <- results %>% left_join(object$chemprop %>% select(cas_number,
                                                              PubChem_CID,
                                                              DTXSID_DTX,
                                                              PREFERRED_NAME,
                                                              CASRN,
                                                              SMILES,
                                                              QSAR_READY_SMILES,
                                                              INCHIKEY,
                                                              MOLECULAR_FORMULA,
                                                              AVERAGE_MASS,
                                                              MONOISOTOPIC_MASS,
                                                              LOG_S,
                                                              LOG_S_AD,
                                                              LOG_S_COMMENT,
                                                              EXCLUDE,
                                                              REMARKS),
                                                  by = "cas_number"
  )

  # reorder list

  results <- results %>% select(EXCLUDE, REMARKS, cas_number, cas, PubChem_CID,
                                dtxsid_ecotox, DTXSID_DTX,
                                PREFERRED_NAME, test_location,
                                reference_number,
                                conc1_type, conc1_mean_op, conc1_mean, conc1_min_op, conc1_min, conc1_max_op, conc1_max,
                                conc1_unit, conc1_comments, concentration_mean, concentration_unit, test_id,
                                endpoint, endpoint_comments, effect, effect_comments,
                                measurement, measurement_comments, obs_duration_mean,
                                obs_duration_unit, organism_lifestage, common_name, latin_name, kingdom,
                                phylum_division, subphylum_div, superclass, class, tax_order, family,
                                genus, species, subspecies, variety, ecotox_group, chemical_name,
                                compound_class, CASRN, SMILES, INCHIKEY, QSAR_READY_SMILES, MOLECULAR_FORMULA, AVERAGE_MASS,
                                MONOISOTOPIC_MASS, LOG_S, LOG_S_AD, LOG_S_COMMENT, reference_db, result_id,
                                reference_type, author, title, source, publication_year, include_endpoint, include_species)
  object$results <- results
  return(object)
}


calculate_hours <- function(object = object){
  message("[EcoToxR]:  Recalculating duration data (day to hours).")
  results <- tibble(object$results)

  results_h <- results %>% filter(obs_duration_unit == "h")
  results_d <- results %>% filter(obs_duration_unit == "d")

  results_d <- results_d %>%
      filter(obs_duration_unit == "d") %>%
      mutate(obs_duration_mean = 24 * obs_duration_mean) %>%
      mutate(obs_duration_unit = "h")

  results <- tibble(bind_rows(results_h, results_d))

  object$results <- results

  return(object)

}


convert_water_solubility <- function(object) {

    message("[EcoToxR]:  Convert the solubility to mg/L.")

    object <- object %>%
        rowwise() %>%
        mutate(S_mg_L = if_else(!is.na(AVERAGE_MASS),
                                true = signif(x = 10^LOG_S * 1000 * AVERAGE_MASS, digits = 4),
                                false = NaN)) %>%
        ungroup() %>%
        rename(QSAR_S_AD = "LOG_S_AD", QSAR_S_COMMENT = "LOG_S_COMMENT")

    return(object)

}



calculate_solubility_domain <- function(object = object, input_column_list = c("quantile_value_mg_L",
                                                                               "min_value_mg_L",
                                                                               "max_value_mg_L",
                                                                               "mean_value_mg_L",
                                                                               "geomean_value_mg_L",
                                                                               "median_value_mg_L")){

    message("[EcoToxR]:  Estimating the solubility domain.")

    # add new dummy columns for the calculation
    #
    #
    #
    #
    #

    for (input_column in input_column_list) {

        object <- object %>% ungroup()

        object <- object %>%
            mutate(input_data = as.vector(pull(object[, input_column]))) %>%
            add_column(input_data_ad = NA) %>%
            rowwise() %>%
            mutate(input_data_ad = case_when(input_data <= S_mg_L ~ 3,

                input_data > S_mg_L & input_data <= 10^5 * log10(S_mg_L) ~ 2,

                input_data > S_mg_L & input_data > 10^5 * log10(S_mg_L) & input_data <= 10^10 * log10(S_mg_L) ~ 1,

                input_data > 10^10 * log10(S_mg_L) ~ 0
                )
            )


        ad_output_column = case_when(input_column == "quantile_value_mg_L" ~ "quantile_value_S_AD",
                                     input_column == "min_value_mg_L" ~ "min_value_S_AD",
                                     input_column == "max_value_mg_L" ~ "max_value_S_AD",
                                     input_column == "mean_value_mg_L" ~ "mean_value_S_AD",
                                     input_column == "geomean_value_mg_L" ~ "geomean_value_S_AD",
                                     input_column == "median_value_mg_L" ~ "median_value_S_AD")


        object <- object %>% select(., -input_data) %>%
            rename(!!ad_output_column := input_data_ad) %>%
            ungroup()

    }


  #results <- object$results

  # Calculate S in mg/L
  #results <- results %>%
  #    rowwise() %>%
  #    mutate(S_mg_L = 10^LOG_S * AVERAGE_MASS * 1000) %>%
  #    filter(is.na(EXCLUDE)) # OPERA output g/L

  # Domain estimate solubility domain (based on ideas in ChemProp)
  # Case 1: if EC <= Sw -> 3
  # Case 2: if 5*log10 Sw >= EC > Sw -> 2
  # Case 3: if 10*log10 Sw >= EC > 5*log10 Sw -> 1
  # Case 4: if EC > 10*log10 Sw -> 0


  #length_progressbar <- length(s)
  # pb <- progress::progress_bar$new(
  #   format = "[EcoToxR]:  Estimating the solubility domain [:bar] :percent ETA: :eta",
  #   total = length_progressbar, clear = FALSE, width = 80)
  #
  # for(i in s){
  # AD <- paste0(i, "_AD")
  # results[, AD] <- NaN
  # # A little bit redundant, but dynamic not yet implemented


# Opera
#




  # Finally exclude if no meaningful data is available

  # results <- results %>%
  #     rowwise() %>%
  #     mutate(
  #             EXCLUDE = case_when(is.na(EXCLUDE) & is.na(concentration_mean) ~ 1)
  #
  #     )

# this is only for debugging
  # # If concentration_mean is <= Sw, the
  # for(j in 1:nrow(results)){
  #  #print(paste0(j, " of ", nrow(results), " in ", i))
  #
  #
  #
  #   if(is.na(results$EXCLUDE[j]) & is.na(results[j, ..i])){
  #     results[j, AD] <- NA
  #
  #   }
  #
  #   else if(is.na(results$EXCLUDE[j]) & is.na(results$concentration_mean[j])){
  #     results[j, AD] <- NA
  #     results$EXCLUDE[j] <- 1
  #   }
  #
  #   else if(is.na(results$EXCLUDE[j]) & results$concentration_mean[j] <= results[j, ..i]){
  #     results[j, AD] <- 3
  #
  #     }
  #
  #   else if(results$concentration_mean[j] > results[j, ..i] & results$concentration_mean[j] <= 10^5 * log10(results[j, ..i])){
  #     results[j, AD] <- 2
  #
  #     }
  #
  #   else if (results$concentration_mean[j] > 10^5 * log10(results[j, ..i]) & results$concentration_mean[j] <= 10^10 * log10(results[j, ..i])){
  #     results[j, AD] <- 1
  #
  #     }
  #
  #   else if (results$concentration_mean[j] > 10^10 * log10(results[j, ..i])){
  #     results[j, AD] <- 0
  #
  #   }
  #   }
  # pb$tick()
  #
  # }
  # #pb$terminate()
  # object$results <- results
  return(object)
}


export_chemical_list <- function(object, project_path){
  message("[EcoToxR]:  Exporting the list of included chemicals to the project folder.")
  message("[EcoToxR]:  Check the file for exclusion of chemicals.")
  message("[EcoToxR]:  Impute missing values for physical-chemical properties.")
  message(paste0("[EcoToxR]:  Edit the file ", tolower(object$parameters$ecotoxgroup), "_chemical_list.csv and re-run the workflow."))

  chemical_list <- object$results_filtered %>%
      select(cas_number, cas, chemical_name) %>%
      unique()

  chemical_list <- object$results_filtered %>%
      select(cas_number, cas, chemical_name) %>%
      unique() %>%
      left_join(object$chemprop %>% select(cas_number, dtxsid_ecotox, PubChem_CID, FOUND_BY, DTXSID_DTX, PREFERRED_NAME, CASRN, INCHIKEY,
                                           IUPAC_NAME, SMILES, INCHI_STRING, MOLECULAR_FORMULA, AVERAGE_MASS,
                                           MONOISOTOPIC_MASS, QSAR_READY_SMILES, QC_LEVEL, LOG_S, LOG_S_AD,
                                           LOG_S_COMMENT, EXCLUDE, REMARKS),
                by = "cas_number") %>%
      arrange(FOUND_BY)

  write_csv(x = chemical_list,
            file = suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")))),
            col_names = TRUE)
}

export_exclude_list <- function(object, project_path){
  message("[EcoToxR]:  Exporting the final list to the project folder.")
  message("[EcoToxR]:  Check the file for exclusion of records.")
  message(paste0("[EcoToxR]:  Edit the file ", tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv and re-run the workflow."))
  exclude <- data.table(matrix(nrow = nrow(object$results_filtered), ncol = 2))
  colnames(exclude) <- c("exclude", "exclusion_comment")
  exclude_list <- data.table(exclude,object$results_filtered)
  object$results_filtered_exclude_list <- data.table(exclude_list[, c("exclude", "exclusion_comment", "result_id", "cas_number", "cas", "chemical_name",
                                                                       colnames(exclude_list)[grep(pattern = "conc1", colnames(exclude_list))],
                                                                       "species_number", "concentration_mean", "concentration_unit",
                                                                       "latin_name", "author", "title", "source", "publication_year"),
                                                                    with = FALSE])
  object$results_filtered_exclude_list <- data.table(object$results_filtered_exclude_list[order(chemical_name)])
  write_csv(data.table(unique(object$results_filtered_exclude_list$cas)), suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list_cas.csv")))), na = "NA")
  write_csv(object$results_filtered_exclude_list, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv")))), na = "NA")
  return(object)
}

remove_excluded_chemicals <- function(object, project_path){

  chemical_list <- read_csv(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")),
                            na = c("NA", "", NA, NaN, "N/A"), show_col_types = FALSE)

  suppressWarnings(chemical_list <- format_chemical_properties(chemical_list))
  object$chemical_list <- chemical_list


  if(!all(unique(chemical_list$EXCLUDE) == 1, na.rm = TRUE)) {
    stop("Only NA or 1 is allowed in the exclusion filter. Check the exclusion list and re-run the workflow.")
  }
  inclusion_list <- chemical_list %>%
      filter(is.na(EXCLUDE)) %>%
      pull(cas_number)

  exclusion_list <- chemical_list %>%
      filter(EXCLUDE == 1) %>%
      pull(cas_number)

  object$results_excluded_by_chemical <- object$results_filtered %>% group_by(cas_number) %>% filter(cas_number %in% exclusion_list)

  object$results_filtered <- object$results_filtered  %>% group_by(cas_number) %>% filter(cas_number %in% inclusion_list)

  write_csv(object$results_filtered, file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_results_included_by_chemical.csv")), na = "NA")

  if(object$results_excluded_by_chemical %>% nrow() > 0){
    write_csv(object$results_excluded, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_results_excluded_by_chemical.csv")))), na = "NA")
  }
  return(object)
}

query_pubchem <- function(object = object) {

    # Split the list in entries with smiles and w/o smiles
    chemical_list_with_SMILES <- object %>% filter(!is.na(SMILES))
    chemical_list_no_SMILES <- object %>% filter(is.na(SMILES))

    # get information from pubchem to fill gaps of DTXSID query
    #
    pubchem <- tibble(
        "cas_number" = integer(),
        "cas" = character(),
        "FOUND_BY" = character(),
        "PREFERRED_NAME" = character(),
        "CASRN" = character(),
        "INCHIKEY" = character(),
        "IUPAC_NAME" = character(),
        "SMILES" = character(),
        "INCHI_STRING" = character(),
        "MOLECULAR_FORMULA" = character(),
        "AVERAGE_MASS" = numeric(),
        "MONOISOTOPIC_MASS" = numeric()
    )


    length_progressbar <- nrow(chemical_list_no_SMILES)
    pb <- progress::progress_bar$new(
        format = "[EcoToxR]:  Retrival of PubChem data [:bar] :percent ETA: :eta",
        total = length_progressbar, clear = FALSE, width = 80)

    for (i in 1:nrow(chemical_list_no_SMILES)){

        pb$tick()


        # lookup for CID based on CASRN
        casrn <- chemical_list_no_SMILES[i, "cas"][[1]]
        pccid <- get_cid(casrn)
        cas_number <- chemical_list_no_SMILES[i, "cas_number"][[1]]
        pccid <- pccid %>% mutate(across(cid, as.integer)) %>% mutate(cas_number = cas_number)


        if (is.na(pccid$cid[[1]])) {
            pubchem_new_row <- tibble(
                "cas_number" = integer(),
                "cas" = character(),
                "FOUND_BY" = character(),
                "PREFERRED_NAME" = character(),
                "CASRN" = character(),
                "INCHIKEY" = character(),
                "IUPAC_NAME" = character(),
                "SMILES" = character(),
                "INCHI_STRING" = character(),
                "MOLECULAR_FORMULA" = character(),
                "AVERAGE_MASS" = numeric(),
                "MONOISOTOPIC_MASS" = numeric()
            )

            cas_numb <- cas_number # rename to avoid confusion

            pubchem_new_row <- pubchem_new_row %>% add_row() %>% mutate(cas_number = cas_numb, cas = casrn, FOUND_BY = "No data retrieved from PubChem")

            pubchem <- bind_rows(pubchem, pubchem_new_row)

            next()

        } else if (nrow(pccid) > 1) {
            # Us the first entry in PubChem only
            pccid <- pccid %>% arrange(cid)
            pccid <- pccid %>% slice_min(cid, n = 1)

        } else {
            pc_props <- tibble(pc_prop(pccid$cid, properties = c("Title",
                                                                 "InChIKey",
                                                                 "IUPACName",
                                                                 "CanonicalSMILES",
                                                                 "InChI",
                                                                 "MolecularFormula",
                                                                 "MolecularWeight",
                                                                 "MonoisotopicMass")))

            if (is.na(pc_props$CanonicalSMILES)) {

                pubchem_new_row <- tibble(
                    "cas_number" = integer(),
                    "cas" = character(),
                    "FOUND_BY" = character(),
                    "PREFERRED_NAME" = character(),
                    "CASRN" = character(),
                    "INCHIKEY" = character(),
                    "IUPAC_NAME" = character(),
                    "SMILES" = character(),
                    "INCHI_STRING" = character(),
                    "MOLECULAR_FORMULA" = character(),
                    "AVERAGE_MASS" = numeric(),
                    "MONOISOTOPIC_MASS" = numeric()
                )

                cas_numb <- cas_number # rename to avoid confusion

                pubchem_new_row <- pubchem_new_row %>% add_row() %>% mutate(cas_number = cas_numb, cas = casrn, FOUND_BY = "No data retrieved from PubChem")

                pubchem <- bind_rows(pubchem, pubchem_new_row)

                next()
            }


            # Postprocess the retrieved data
            pc_props <- pccid %>% left_join(pc_props, by = c("cid" = "CID"))

            if (!is_empty(which(names(pc_props) %like% "IUPACName"))) {

                pc_props <- pc_props %>% rename(c("cas" = "query", "cid" = "cid", "cas_number" = "cas_number",
                                                  "PREFERRED_NAME" = "Title",
                                                  "MOLECULAR_FORMULA" = "MolecularFormula",
                                                  "AVERAGE_MASS" = "MolecularWeight", "SMILES" = "CanonicalSMILES",
                                                  "INCHI_STRING"  = "InChI", "INCHIKEY" = "InChIKey",
                                                  "IUPAC_NAME" = "IUPACName", "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

                pc_props <- pc_props %>% mutate("CASRN" = cas)

                pc_props <- pc_props %>% select(cas_number, cas, PREFERRED_NAME, CASRN,
                                                INCHIKEY, IUPAC_NAME, SMILES, INCHI_STRING, MOLECULAR_FORMULA,
                                                AVERAGE_MASS, MONOISOTOPIC_MASS)

            } else {

                pc_props <- pc_props %>% rename(c("cas" = "query", "cid" = "cid", "cas_number" = "cas_number",
                                                  "PREFERRED_NAME" = "Title",
                                                  "MOLECULAR_FORMULA" = "MolecularFormula",
                                                  "AVERAGE_MASS" = "MolecularWeight", "SMILES" = "CanonicalSMILES",
                                                  "INCHI_STRING"  = "InChI", "INCHIKEY" = "InChIKey",
                                                  "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

                pc_props <- pc_props %>% mutate(IUPAC_NAME = NA)
                pc_props <- pc_props %>% mutate(CASRN = cas)

                pc_props <- pc_props %>% select(cas_number, cas, PREFERRED_NAME, CASRN,
                                                INCHIKEY, IUPAC_NAME, SMILES, INCHI_STRING, MOLECULAR_FORMULA,
                                                AVERAGE_MASS, MONOISOTOPIC_MASS)
            }

            pc_props <- pc_props %>% add_column(FOUND_BY = "PubChem", .before = 3)


        }

        pubchem <- pubchem %>% rbind(pc_props)

    }

    # Output of webchem is character, needs to be fixed here.
    pubchem <- pubchem %>% mutate(across(cas_number, as.integer)) %>% mutate(across(AVERAGE_MASS:MONOISOTOPIC_MASS, as.double))

    chemicals_update <- chemical_list_no_SMILES %>% inner_join(pubchem %>% select(cas_number))

    chemicals_update <- chemicals_update %>%
        mutate(FOUND_BY = pubchem$FOUND_BY, PREFERRED_NAME = pubchem$PREFERRED_NAME, CASRN = pubchem$CASRN,
               INCHIKEY = pubchem$INCHIKEY, IUPAC_NAME = pubchem$IUPAC_NAME, SMILES = pubchem$SMILES,
               INCHI_STRING = pubchem$INCHI_STRING, MOLECULAR_FORMULA = pubchem$MOLECULAR_FORMULA,
               AVERAGE_MASS = pubchem$AVERAGE_MASS, MONOISOTOPIC_MASS = pubchem$MONOISOTOPIC_MASS
        )


    chemical_list_no_SMILES <- chemical_list_no_SMILES %>% filter(cas_number != chemicals_update$cas_number)


    # recombine lists
    chemical_list <- bind_rows(chemical_list_with_SMILES, chemicals_update)


    # add comments and exclude those entries with remaining gaps
    chemical_list <- chemical_list %>% mutate(EXCLUDE = if_else(condition = FOUND_BY == "No data retrieved from PubChem",
                                                                true = 1,
                                                                false = EXCLUDE),
                                              REMARKS = if_else(condition = FOUND_BY == "No data retrieved from PubChem",
                                                                true = "no data",
                                                                false = REMARKS))


}

update_chemical_list <- function(object = object, project_path = project_path, database_path = database_path){

  # some kind dirty solution...
  #object$chemprop <- object$chemprop[chemprop, on = c("cas_number"), AVERAGE_MASS := i.AVERAGE_MASS]
  #object$chemprop <- object$chemprop %>% mutate_all(na_if,"")
  #

  chemical_list_extended <- read_csv(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")),
                              na = c("NA", "", NA, NaN, "N/A"), show_col_types = FALSE)

  chemprop <- data.table(object$chemprop)
  chemical_list <- data.table(object$chemical_list)
  chemical_list <- format_chemical_properties(chemical_list)
  chemprop <- data.table(merge(chemprop, chemical_list, by = 'cas_number', all = TRUE, suffixes = c("", ".update")))
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
chemprop[, (col_names2) := NULL]
export_chemical_properties(object = chemprop, database_path = database_path, project_path = project_path)
object$chemprop <- tibble(chemprop)

return(object)
}

# Handle the chemical lists
export_chemical_properties <- function(object, database_path = database_path,
                                       project_path = project_path){

  object <- format_chemical_properties(object)
  object <- tibble(object)

  write_csv(x = object, file = suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))), na = "NA")
  write_csv(x = object, file = suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))), na = "NA")

}

create_chemical_properties <- function(database_path){
  if (file.exists(file.path(database_path, "chemical_properties.csv"))){
    message("[EcoToxR]:  Reading chemical properties (this is a custom file).")
    object <- read_csv(file = file.path(database_path, "chemical_properties.csv"), na = c("NA", "", NA, NaN), show_col_types = FALSE)
    object <- format_chemical_properties(object)

  } else {

    message("[EcoToxR]:  A custom file for the storage of chemical properties (mol weight) is compiled.")
    message("[EcoToxR]:  The file is stored in the current database folder and will be copied to the project folder for backup.")

    object <- tibble(
        "cas_number" = integer(),
        "cas" = character(),
        "chemical_name" = character(),
        "dtxsid_ecotox" = character(),
        "PubChem_CID" = integer(),
        "FOUND_BY" = character(),
        "DTXSID_DTX" = character(),
        "PREFERRED_NAME" = character(),
        "CASRN" = character(),
        "INCHIKEY"  = character(),
        "IUPAC_NAME" = character(),
        "SMILES" = character(),
        "INCHI_STRING" = character(),
        "QSAR_READY_SMILES" = character(),
        "MOLECULAR_FORMULA" = character(),
        "AVERAGE_MASS" = numeric(),
        "MONOISOTOPIC_MASS" = numeric(),
        "QC_LEVEL" = integer(),
        "LOG_S" = numeric(),
        "LOG_S_AD" = integer(),
        "LOG_S_COMMENT" = character(),
        "EXCLUDE" = integer(),
        "REMARKS" = character()
    )

    # object <- format_chemical_properties(object)
    suppressWarnings(
      write_csv(x = object, file = suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))),
                na = "NA",
                col_names = TRUE)
    )

  }
  return(object)
}

format_chemical_properties <- function(object){
      suppressWarnings(object$cas_number <- as.integer(object$cas_number))
      object$cas <- as.character(object$cas)
      object$chemical_name <- as.character(object$chemical_name)
      object$dtxsid_ecotox <- as.character(object$dtxsid_ecotox)
      object$FOUND_BY <-  as.character(object$FOUND_BY)
      object$DTXSID_DTX <- as.character(object$DTXSID_DTX)
      object$PubChem_CID <- as.integer(object$PubChem_CID)
      object$PREFERRED_NAME <- as.character(object$PREFERRED_NAME)
      object$CASRN <- as.character(object$CASRN)
      object$INCHIKEY <- as.character(object$INCHIKEY)
      object$IUPAC_NAME <- as.character(object$IUPAC_NAME)
      object$SMILES <- as.character(object$SMILES)
      object$INCHI_STRING <- as.character(object$INCHI_STRING)
      object$MOLECULAR_FORMULA <- as.character(object$MOLECULAR_FORMULA)
      suppressWarnings(object$AVERAGE_MASS <- as.numeric(object$AVERAGE_MASS))
      suppressWarnings(object$MONOISOTOPIC_MASS <- as.numeric(object$MONOISOTOPIC_MASS))
      object$QSAR_READY_SMILES <- as.character(object$QSAR_READY_SMILES)
      object$QC_LEVEL <- as.integer(object$QC_LEVEL)
      object$LOG_S <- as.numeric(object$LOG_S)
      object$LOG_S_AD <- as.integer(object$LOG_S_AD)
      object$LOG_S_COMMENT <- as.character(object$LOG_S_COMMENT)
      object$EXCLUDE <- as.numeric(object$EXCLUDE)
      object$REMARKS <- as.character(object$REMARKS)
  return(object)
}


export_mol_units <- function(object, project_path = project$project_path) {

  ecotoxgroup <- object$parameters$ecotoxgroup
  file_name <- suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))))

  message("[EcoToxR]:  Exporting the list of chemicals for mol/L to mg/L conversion to the project folder.")
  message(paste0("[EcoToxR]:  Edit the file ", tolower(ecotoxgroup), "_mol_weight.csv and re-run the workflow."))
  message("[EcoToxR]:  Add only missing data to column 'AVERAGE_MASS'.")

  chemprop <- object$results_filtered %>%
      filter(conc1_unit %like% "mol/L") %>%
      select(cas_number, cas, chemical_name, dtxsid_ecotox) %>%
      add_column("AVERAGE_MASS" = NA)

  if(nrow(object$chemprop) > 0) {

      chemprop <- chemprop %>%
          left_join(object$chemprop %>%
                        select(1, 4:ncol(object$chemprop)), by = "cas_number")
      }

  chemprop$AVERAGE_MASS.x <- NULL

  col_chemprop <- colnames(chemprop)
  col_num <- which(col_chemprop == "AVERAGE_MASS.y")
  colnames(chemprop)[col_num] <- "AVERAGE_MASS"


  chemprop <- chemprop[order(chemprop$chemical_name), ]

  chemprop <- chemprop %>%
      group_by(cas_number) %>%
      unique()


  #chemprop <- chemprop %>%
  #    group_by(FOUND_BY) %>%
  #    arrange(desc(FOUND_BY))

  write_csv(x = chemprop, file = file_name, col_names = TRUE)

  return(object)
}

update_mol_units <- function(object, database_path = project$database_path, project_path = project$project_path) {
  message("[EcoToxR]:  Recalulating record unit mol/L to mg/L.")
  ecotoxgroup <- object$parameters$ecotoxgroup
  file_name <- file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))

  chemprop <- read_csv(file = file_name, col_names = TRUE, na = c("NA", "", NA, NaN, "N/A", "n/a"), show_col_types = FALSE)

  # chemprop <- format_chemical_properties(chemprop)

  # subset

  mol_records <- object$results_filtered %>%
      filter(concentration_unit %like% "mol/L") %>%
      left_join(chemprop %>% select(cas_number, AVERAGE_MASS), by = "cas_number") %>%
      rowwise() %>%
      mutate(concentration_mean = concentration_mean * 1000 * AVERAGE_MASS) %>%
      mutate(concentration_unit = "mg/L") %>%
      select(-AVERAGE_MASS)

  object$results_filtered <- rbind(object$results_filtered %>% filter(concentration_unit == "mg/L"), mol_records)

  message("[EcoToxR]: The data finally contains the following units ")
  print(unique(object$results_filtered$concentration_unit))

  #message("[EcoToxR]: Update the chemical properties.")

  object$chemprop <- tibble(data.table(object$chemprop)[data.table(chemprop), on = c("cas_number"), AVERAGE_MASS := i.AVERAGE_MASS])

  object$chemprop <- object$chemprop %>% mutate_all(na_if, "")

  write_csv(x = object$chemprop, file = suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))), col_names = TRUE)

  return(object)
}

convert_units <- function(object, sample_size = NA) {

  if(!is.na(sample_size)) {
    object <- object %>% sample_n(size = sample_size)
  }

  message("[EcoToxR]:  Impute missing concentration mean by averaging minimum and maximum values")

  object <- object %>% mutate(concentration_mean = conc1_mean, concentration_unit = conc1_unit)

  object <- bind_rows(object %>% filter(!is.na(concentration_mean)),
                  object %>% filter(is.na(concentration_mean)) %>%
                      rowwise() %>%
                      mutate(concentration_mean = if_else(condition = is.na(concentration_mean),
                                                          true = mean(c(conc1_min, conc1_max)),
                                                          false = concentration_mean))
  )

#
#
#
#
#   length_progressbar <- length(object$concentration_mean[is.na(object$concentration_mean)])
#   if(length_progressbar > 0){
#     pb <- progress::progress_bar$new(
#       format = "[EcoToxR]:  Update missing concentration mean with mean from min and max (:spin) [:bar] :percent ETA: :eta",
#       total = length_progressbar, clear = FALSE, width = 80)
#
#     for(i in 1:length(object$concentration_mean[is.na(object$concentration_mean)])) {
#       pb$tick()
#       object[which(object[, is.na(concentration_unit)])][i]$concentration_mean <-
#         mean(object[which(object[, is.na(concentration_unit)])][i]$conc1_min, object[which(object[, is.na(concentration_unit)])][i]$conc1_max)
#
#
#

  # convert similar unit to SI conform units (e.g. ppm to mg/L)

  # ng related
  message("[EcoToxR]:  Converting units to human readible format (i.e. mg/L).")
  message("[EcoToxR]:  The following units will be converted:")
  print(unique(object$concentration_unit))

  # this assumes that only some records are reported in log(LC50)
  log <- "^.log.[A-Z]*[0-9]*?$"

  if (nrow(object %>% filter(endpoint %like% log)) > 0) {

      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = endpoint %like% log,
                                                                             true = 10^concentration_mean,
                                                                             false= concentration_mean))
      object <- object %>% rowwise() %>% mutate(endpoint = if_else(condition = concentration_unit %like% log,
                                                                   true = "LC50",
                                                                   false = endpoint))

  }

  # pg/L
  pg <- "^pg/(l|L)$"

  if (nrow(object %>% filter(concentration_unit %like% pg)) > 0) {

      message("[EcoToxR]:  Converting pg/L like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% pg,
                                                                             true = concentration_mean / 1e+09,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% pg,
                                                                             true = "mg/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping pg/L like units")}


  #ng
  ng <- "^(|A(I|E|i|e) )ng/(l|L|dm3)$|^ppt$|^(|A(I|E|i|e) )pg/(ml|mL)$"

  if (nrow(object %>% filter(concentration_unit %like% ng)) > 0) {

      message("[EcoToxR]:  Converting ng/L like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% ng,
                                                                             true = concentration_mean / 1e+06,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% ng,
                                                                             true ="mg/L",
                                                                             false =concentration_unit))

  } else {message("[EcoToxR]:  Skipping ng/L like units")}

  #µg related
  ug <- "(^(|(A|a)(I|E|i|e) )ug/(l|L|dm3)$)|^mg/m3$|^ppb$|(^pg/u(l|L)$)|(^(|(A|a)(I|E|i|e) )ng/m(l|L)$)|^pg/u(l|L)"

  if (nrow(object %>% filter(concentration_unit %like% ug)) > 0) {

      message("[EcoToxR]:  Converting ug/L like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% ug,
                                                                             true = concentration_mean / 1e+03,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% ug,
                                                                             true = "mg/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping ug/L like units")}


  # g/L
  g <- "(^(|(A|a)(e|i|E|I) )(g/(l|L)$))|(^mg/(ml|mL)$|(^g/dm3$))|(^ug/(ul|uL)$)"

  if (nrow(object %>% filter(concentration_unit %like% g)) > 0) {

      message("[EcoToxR]:  Converting g/L like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% g,
                                                                             true = concentration_mean * 1e+03,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% g,
                                                                             true = "mg/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping g/L like units")}


  # mg/L
  message("[EcoToxR]:  Converting mg/L like units")
  mg <- "(^(|((A|a)(I|E|i|e)) )mg/(l|L|dm3)$)|^ppm$|(^(|((A|a)(I|E|i|e)) )ug/m(l|L|m3)$)|^g/m3$"
  object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% mg,
                                                                         true = "mg/L",
                                                                         false = concentration_unit))

  #pmol related
  pmol <- "(^(|(A|a)(I|E|i|e) )pmol/(l|L)$)"

  if (nrow(object %>% filter(concentration_unit %like% pmol)) > 0) {

      message("[EcoToxR]:  Converting pmol like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% pmol,
                                                                             true = concentration_mean / 1e+12,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% pmol,
                                                                            true = "mol/L",
                                                                            false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping pmol/L like units")}

  #nmol related
  nmol <- "(^(|(A|a)(I|E|i|e) )nmol/(l|L)$)"

  if (nrow(object %>% filter(concentration_unit %like% nmol)) > 0) {

      message("[EcoToxR]:  Converting nmol like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% nmol,
                                                                             true = concentration_mean / 1e+09,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% nmol,
                                                                             true = "mol/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping nmol/L like units")}

  #umol related
  umol <- "(^(|(A|a)(I|E|i|e) )umol/(dm3|l|L)$)|^nmol/ml$|^mmol/m3$"

  if (nrow(object %>% filter(concentration_unit %like% umol)) > 0) {

      message("[EcoToxR]:  Converting pmol like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% umol,
                                                                             true = concentration_mean / 1e+06,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% umol,
                                                                             true = "mol/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping umol/L like units")}

  #mmol related
  mmol <- "^(|(A|a)(I|E|i|e) )mmol/(dm3|l|L)$"

  if (nrow(object %>% filter(concentration_unit %like% mmol)) > 0) {

      message("[EcoToxR]:  Converting mmol like units")
      object <- object %>% rowwise() %>% mutate(concentration_mean = if_else(condition = concentration_unit %like% mmol,
                                                                             true = concentration_mean / 1e+03,
                                                                             false = concentration_mean))
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% mmol,
                                                                             true = "mol/L",
                                                                             false = concentration_unit))

  } else {message("[EcoToxR]:  Skipping mmol/L like units")}


  #mol related
  mol <- "^(|(A|a)(I|E|i|e) )mol/(dm3)$"

  if (nrow(object %>% filter(concentration_unit %like% mol)) > 0) {
      message("[EcoToxR]:  Converting mol like units")
      object <- object %>% rowwise() %>% mutate(concentration_unit = if_else(condition = concentration_unit %like% mol,
                                                                             true = "mol/L",
                                                                             false = concentration_unit))
  } else {message("[EcoToxR]:  Skipping mol/L like units")}


  message("[EcoToxR]:  The following units are left in the data:")
  print(unique(object$concentration_unit))


  unit_test <- "(mg/L|mol/L)"

  if (any(!unique(object$concentration_unit %like% unit_test))) {
      message("[EcoToxR]")
      message("[EcoToxR]:  ##########################################################")
      message("[EcoToxR]:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      message("[EcoToxR]:  The remaining units contain other than mg/L and mol/L.   !")
      message("[EcoToxR]:  Contact the developer to update the unit conversion.     !")
      message("[EcoToxR]:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      message("[EcoToxR]:  ##########################################################")
      message("[EcoToxR]")
  }

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

save_project <- function(object = object, save_project_steps = save_project_steps){
  state <- object$object$state
  project_path <- object$project_path
  ecotoxgroup <- object$object$parameters$ecotoxgroup
  if(isTRUE(save_project_steps)){
    project <- object
    message("[EcoToxR]:  Saving the project state.")
    save(project,
         file = file.path(project_path, paste0(tolower(ecotoxgroup), "_state", state, ".RData")),
         compress = TRUE)
  }
}

get_git <- function(){


  #chemical_list <- read_csv("path_to_project/algae_chemical_list.csv")

  length_progressbar <- nrow(chemical_list)
  pb <- progress::progress_bar$new(
    format = "[EcoToxR]:  Retrival of PubChem data [:bar] :percent ETA: :eta",
    total = length_progressbar, clear = FALSE, width = 80)

  for (i in 1:nrow(chemical_list)){

    pb$tick()

    if (is.na(chemical_list$PubChem_CID[i])){
      # lookup for CID based on CASRN
      casrn <- chemical_list[i, "CASRN"][[1]]
      pccid <- get_cid(casrn)
      cas_number <- chemical_list[i, "cas_number"][[1]]
      pccid <- pccid %>% mutate(across(CID, as.integer)) %>% mutate(cas_number = cas_number)


      if (is.na(pccid$CID[[1]])) {
        next()

      } else if (nrow(pccid) > 1) {
        # Us the first entry in PubChem only
        pccid <- pccid %>% arrange(PubChem_CID)
        pccid <- pccid %>% slice_min(PubChem_CID, n = 1)

      }


      chemical_list[i, "CID"] <- pccid[[2]]

      # write_csv(chemical_list, "c:/TEMP/EcoToxDB/test/algae_chemical_list.csv")

    } else {
      next()
    }

  }


  write_csv(chemical_list, "path_to_project/algae_chemical_list.csv")

}
