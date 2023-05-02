# helper functions (functions without real allocation)

remove_excluded_chemicals <- function(object, project_path){

  chemical_list <- readr::read_csv(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")),
                            na = c("NA", "", NA, NaN, "N/A"), show_col_types = FALSE)

  suppressWarnings(chemical_list <- format_chemical_properties(chemical_list))
  object$chemical_list <- chemical_list


  if(!all(unique(chemical_list$EXCLUDE) == 1, na.rm = TRUE)) {
    stop("Only NA or 1 is allowed in the exclusion filter. Check the exclusion list and re-run the workflow.")
  }
  inclusion_list <- chemical_list %>%
    dplyr::filter(is.na(EXCLUDE)) %>%
    dplyr::pull(cas_number) %>%
    as.double()

  exclusion_list <- chemical_list %>%
    dplyr::filter(EXCLUDE == 1) %>%
    dplyr::pull(cas_number) %>%
    as.double()

  object$results_excluded_by_chemical <-
    object$results_filtered %>%
    dplyr::group_by(cas_number) %>%
    dplyr::filter(cas_number %in% exclusion_list)

  object$results_filtered <-
    object$results_filtered %>%
    dplyr::group_by(cas_number) %>%
    dplyr::filter(cas_number %in% inclusion_list)

  readr::write_csv(object$results_filtered, file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_results_included_by_chemical.csv")), na = "NA")

  if(object$results_excluded_by_chemical %>% nrow() > 0){
    readr::write_csv(object$results_excluded, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_results_excluded_by_chemical.csv")))), na = "NA")
  }
  return(object)
}

# no yet implemented

query_pubchem <- function(object = object) {

    # Split the list in entries with smiles and w/o smiles
    chemical_list_with_SMILES <-
      object %>% dplyr::filter(!is.na(SMILES))

    chemical_list_no_SMILES <-
      object %>% dplyr::filter(is.na(SMILES))

    # get information from pubchem to fill gaps of DTXSID query
    #
    pubchem <- tibble(
        "cas_number" = double(),
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

        pccid <- webchem::get_cid(casrn)

        cas_number <-
          chemical_list_no_SMILES[i, "cas_number"][[1]] %>% as.double()

        pccid <- pccid %>%
          dplyr::mutate(dplyr::across(cid, as.double)) %>%
          dplyr::mutate(cas_number = cas_number)


        if (is.na(pccid$cid[[1]])) {
            pubchem_new_row <- tibble::tibble(
                "cas_number" = double(),
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

            pubchem_new_row <-
              pubchem_new_row %>%
              dplyr::add_row() %>%
              dplyr::mutate(cas_number = cas_numb,
                            cas = casrn,
                            FOUND_BY = "No data retrieved from PubChem")

            pubchem <- dplyr::bind_rows(pubchem, pubchem_new_row)

            next()

        } else if (nrow(pccid) > 1) {
            # Us the first entry in PubChem only
            pccid <- pccid %>% dplyr::arrange(cid)
            pccid <- pccid %>% dplyr::slice_min(cid, n = 1)

        } else {
            pc_props <-
              tibble::tibble(pc_prop(pccid$cid,
                                     properties = c("Title",
                                                    "InChIKey",
                                                    "IUPACName",
                                                    "CanonicalSMILES",
                                                    "InChI",
                                                    "MolecularFormula",
                                                    "MolecularWeight",
                                                    "MonoisotopicMass")))

            if (is.na(pc_props$CanonicalSMILES)) {

                pubchem_new_row <- tibble(
                    "cas_number" = double(),
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

                pubchem_new_row <-
                  pubchem_new_row %>%
                  dplyr::add_row() %>%
                  dplyr::mutate(cas_number = cas_numb,
                                cas = casrn,
                                FOUND_BY = "No data retrieved from PubChem")

                pubchem <- dplyr::bind_rows(pubchem, pubchem_new_row)

                next()
            }


            # Postprocess the retrieved data
            pc_props <-
              pccid %>% dplyr::left_join(pc_props, by = c("cid" = "CID"))

            if (!purrr::is_empty(which(names(pc_props) %like% "IUPACName"))) {

                pc_props <-
                  pc_props %>%
                  dplyr::rename(c("cas" = "query",
                                  "cid" = "cid",
                                  "cas_number" = "cas_number",
                                  "PREFERRED_NAME" = "Title",
                                  "MOLECULAR_FORMULA" = "MolecularFormula",
                                  "AVERAGE_MASS" = "MolecularWeight",
                                  "SMILES" = "CanonicalSMILES",
                                  "INCHI_STRING"  = "InChI",
                                  "INCHIKEY" = "InChIKey",
                                  "IUPAC_NAME" = "IUPACName",
                                  "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

                pc_props <- pc_props %>% dplyr::mutate("CASRN" = cas)

                pc_props <-
                  pc_props %>%
                  dplyr::select(cas_number,
                                cas,
                                PREFERRED_NAME,
                                CASRN,
                                INCHIKEY,
                                IUPAC_NAME,
                                SMILES,
                                INCHI_STRING,
                                MOLECULAR_FORMULA,
                                AVERAGE_MASS,
                                MONOISOTOPIC_MASS)

            } else {

                pc_props <-
                  pc_props %>%
                  dplyr::rename(c("cas" = "query",
                                  "cid" = "cid",
                                  "cas_number" = "cas_number",
                                  "PREFERRED_NAME" = "Title",
                                  "MOLECULAR_FORMULA" = "MolecularFormula",
                                  "AVERAGE_MASS" = "MolecularWeight",
                                  "SMILES" = "CanonicalSMILES",
                                  "INCHI_STRING"  = "InChI",
                                  "INCHIKEY" = "InChIKey",
                                  "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

                pc_props <- pc_props %>% dplyr::mutate(IUPAC_NAME = NA)
                pc_props <- pc_props %>% dplyr::mutate(CASRN = cas)

                pc_props <-
                  pc_props %>% dplyr::select(cas_number,
                                             cas,
                                             PREFERRED_NAME,
                                             CASRN,
                                             INCHIKEY,
                                             IUPAC_NAME,
                                             SMILES,
                                             INCHI_STRING,
                                             MOLECULAR_FORMULA,
                                             AVERAGE_MASS,
                                             MONOISOTOPIC_MASS)
            }

            pc_props <-
              pc_props %>% tibble::add_column(FOUND_BY = "PubChem", .before = 3)


        }

        pubchem <- pubchem %>% dplyr::bind_rows(pc_props)

    }

    # Output of webchem is character, needs to be fixed here.
    pubchem <-
      pubchem %>%
      dplyr::mutate(dplyr::across(cas_number, as.double)) %>%
      dplyr::mutate(across(AVERAGE_MASS:MONOISOTOPIC_MASS, as.double))

    chemicals_update <-
      chemical_list_no_SMILES %>%
      dplyr::inner_join(pubchem %>% dplyr::select(cas_number))

    chemicals_update <- chemicals_update %>%
      dplyr::mutate(FOUND_BY = pubchem$FOUND_BY,
                    PREFERRED_NAME = pubchem$PREFERRED_NAME,
                    CASRN = pubchem$CASRN,
                    INCHIKEY = pubchem$INCHIKEY,
                    IUPAC_NAME = pubchem$IUPAC_NAME,
                    SMILES = pubchem$SMILES,
                    INCHI_STRING = pubchem$INCHI_STRING,
                    MOLECULAR_FORMULA = pubchem$MOLECULAR_FORMULA,
                    AVERAGE_MASS = pubchem$AVERAGE_MASS,
                    MONOISOTOPIC_MASS = pubchem$MONOISOTOPIC_MASS
        )


    chemical_list_no_SMILES <-
      chemical_list_no_SMILES %>%
      dplyr::filter(cas_number != chemicals_update$cas_number)


    # recombine lists
    chemical_list <-
      dplyr::bind_rows(chemical_list_with_SMILES, chemicals_update)


    # add comments and exclude those entries with remaining gaps
    chemical_list <-
      chemical_list %>%
      dplyr::mutate(EXCLUDE =
                      dplyr::if_else(condition = FOUND_BY == "No data retrieved from PubChem",
                                     true = 1,
                                     false = EXCLUDE),
                    REMARKS =
                      dplyr::if_else(condition = FOUND_BY == "No data retrieved from PubChem",
                                     true = "no data",
                                     false = REMARKS))
}

update_chemical_list <- function(object = object,
                                 project_path = project_path,
                                 database_path = database_path) {

  chemical_list_extended <-
    readr::read_csv(file.path(project_path,
                              paste0(tolower(object$parameters$ecotoxgroup),
                                     "_chemical_list.csv")),
                              na = c("NA", "", NA, NaN, "N/A"),
                    show_col_types = FALSE)

  chemprop <- data.table::data.table(object$chemprop)

  chemical_list <- data.table::data.table(object$chemical_list)

  chemical_list <- format_chemical_properties(chemical_list)

  chemprop <-
    data.table::data.table(merge(chemprop,
                                 chemical_list,
                                 by = 'cas_number',
                                 all = TRUE,
                                 suffixes = c("", ".update")))

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

chemprop <- data.table::data.table(chemprop)

chemprop[, (col_names2) := NULL]

export_chemical_properties(object = chemprop,
                           database_path = database_path,
                           project_path = project_path)

object$chemprop <- tibble::tibble(chemprop)

return(object)
}

# Handle the chemical lists


create_chemical_properties <- function(database_path){
  if (file.exists(file.path(database_path, "chemical_properties.csv"))){
    message("[EcoToxR]:  Reading chemical properties (this is a custom file).")
    object <- readr::read_csv(file = file.path(database_path, "chemical_properties.csv"), na = c("NA", "", NA, NaN), show_col_types = FALSE)
    object <- format_chemical_properties(object)

  } else {

    message("[EcoToxR]:  A custom file for the storage of chemical properties (mol weight) is compiled.")
    message("[EcoToxR]:  The file is stored in the current database folder and will be copied to the project folder for backup.")

    object <- tibble::tibble(
        "cas_number" = double(),
        "cas" = character(),
        "chemical_name" = character(),
        "dtxsid_ecotox" = character(),
        "PubChem_CID" = double(),
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
        "S_AD_index" = numeric(),
        "LOG_S_COMMENT" = character(),
        "EXCLUDE" = integer(),
        "REMARKS" = character()
    )

    # object <- format_chemical_properties(object)
    suppressWarnings(
      readr::write_csv(x = object,
                       file = suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))),
                       na = "NA",
                       col_names = TRUE)
    )

  }
  return(object)
}

format_chemical_properties <- function(object){
      object$cas_number <- as.double(object$cas_number)
      object$cas <- as.character(object$cas)
      object$chemical_name <- as.character(object$chemical_name)
      object$FOUND_BY <-  as.character(object$FOUND_BY)
      object$DTXSID_DTX <- as.character(object$DTXSID_DTX)
      object$PubChem_CID <- as.double(object$PubChem_CID)
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
      object$S_AD_index <- as.numeric(object$S_AD_index)
      object$LOG_S_COMMENT <- as.character(object$LOG_S_COMMENT)
      object$EXCLUDE <- as.numeric(object$EXCLUDE)
      object$REMARKS <- as.character(object$REMARKS)
  return(object)
}


get_cid <- function(){

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
      pccid <- webchem::get_cid(casrn)
      cas_number <- chemical_list[i, "cas_number"][[1]] %>% as.double()
      pccid <-
        pccid %>%
        dplyr::mutate(dplyr::across(CID, as.double)) %>%
        dplyr::mutate(cas_number = cas_number)


      if (is.na(pccid$CID[[1]])) {
        next()

      } else if (nrow(pccid) > 1) {
        # Us the first entry in PubChem only
        pccid <- pccid %>% dplyr::arrange(PubChem_CID)
        pccid <- pccid %>% dplyr::slice_min(PubChem_CID, n = 1)

      }


      chemical_list[i, "CID"] <- pccid[[2]]

      # write_csv(chemical_list, "c:/TEMP/EcoToxDB/test/algae_chemical_list.csv")

    } else {
      next()
    }

  }


  write_csv(chemical_list, "path_to_project/algae_chemical_list.csv")

}
