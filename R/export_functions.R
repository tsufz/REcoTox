
# Export functions
export_mol_units <- function(object, project_path = project$project_path) {

    ecotoxgroup <- object$parameters$ecotoxgroup
    file_name <- suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))))

    message("[EcoToxR]:  Exporting the list of chemicals for mol/L to mg/L conversion to the project folder.")
    message(paste0("[EcoToxR]:  Edit the file ", tolower(ecotoxgroup), "_mol_weight.csv and re-run the workflow."))
    message("[EcoToxR]:  Add only missing data to column 'AVERAGE_MASS'.")

    chemprop <- object$results_filtered %>%
        dplyr::filter(conc1_unit %like% "mol/L") %>%
        dplyr::select(cas_number, cas, chemical_name, dtxsid_ecotox) %>%
        tibble::add_column("AVERAGE_MASS" = NA)

    if(nrow(object$chemprop) > 0) {

        chemprop <- chemprop %>%
            dplyr::left_join(object$chemprop %>%
                                 dplyr::select(1, 4:ncol(object$chemprop)),
                             by = "cas_number")
    }

    chemprop$AVERAGE_MASS.x <- NULL

    col_chemprop <- colnames(chemprop)
    col_num <- which(col_chemprop == "AVERAGE_MASS.y")
    colnames(chemprop)[col_num] <- "AVERAGE_MASS"


    chemprop <- chemprop[order(chemprop$chemical_name), ]

    chemprop <- chemprop %>%
        dplyr::group_by(cas_number) %>%
        unique() %>%
        dplyr::ungroup()

    readr::write_csv(x = chemprop,
                     file = file_name,
                     col_names = TRUE)

    return(object)
}

# Export the chemical list for review
export_chemical_list <- function(object, project_path){
    message("[EcoToxR]:  Exporting the list of included chemicals to the project folder.")
    message("[EcoToxR]:  Check the file for exclusion of chemicals.")
    message("[EcoToxR]:  Impute missing values for physical-chemical properties.")
    message(paste0("[EcoToxR]:  Edit the file ", tolower(object$parameters$ecotoxgroup), "_chemical_list.csv and re-run the workflow."))

    chemical_list <- object$results_filtered %>%
        dplyr::select(cas_number, cas, chemical_name) %>%
        unique()

    chemical_list <- object$results_filtered %>%
        dplyr::select(cas_number, cas, chemical_name) %>%
        unique() %>%
        dplyr::left_join(object$chemprop %>%
                             dplyr::select(cas_number, dtxsid_ecotox, PubChem_CID,
                                           FOUND_BY, DTXSID_DTX, PREFERRED_NAME,
                                           CASRN, INCHIKEY, IUPAC_NAME, SMILES,
                                           INCHI_STRING, MOLECULAR_FORMULA,
                                           AVERAGE_MASS, MONOISOTOPIC_MASS,
                                           QSAR_READY_SMILES, QC_LEVEL, LOG_S,
                                           LOG_S_AD, S_AD_index, LOG_S_COMMENT,
                                           EXCLUDE, REMARKS),
                         by = "cas_number") %>%
        dplyr::arrange(FOUND_BY)

    readr::write_csv(x = chemical_list,
                     file = suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_chemical_list.csv")))),
                     col_names = TRUE)
}

# Export the exclusion list
export_exclude_list <- function(object, project_path){
    message("[EcoToxR]:  Exporting the final list to the project folder.")
    message("[EcoToxR]:  Check the file for exclusion of records.")
    message(paste0("[EcoToxR]:  Edit the file ", tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv and re-run the workflow."))

    exclude <-
        data.table(matrix(nrow = nrow(object$results_filtered), ncol = 2))

    colnames(exclude) <-
        c("exclude", "exclusion_comment")

    exclude_list <-
        data.table::data.table(exclude, object$results_filtered)

    object$results_filtered_exclude_list <-
        data.table::data.table(exclude_list[, c("exclude", "exclusion_comment", "result_id", "cas_number", "cas", "chemical_name",
                                                colnames(exclude_list)[grep(pattern = "conc1", colnames(exclude_list))],
                                                "species_number", "concentration_mean", "concentration_unit",
                                                "latin_name", "author", "title", "source", "publication_year"),
                                            with = FALSE])

    object$results_filtered_exclude_list <-
        data.table::data.table(object$results_filtered_exclude_list[order(chemical_name)])

    readr::write_csv(data.table(unique(object$results_filtered_exclude_list$cas)), suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list_cas.csv")))), na = "NA")

    readr::write_csv(object$results_filtered_exclude_list, suppressWarnings(normalizePath(file.path(project_path, paste0(tolower(object$parameters$ecotoxgroup), "_mortality_filtered_exclude_list.csv")))), na = "NA")

    return(object)
}

# Export the chemical properties in the database folder
export_chemical_properties <- function(object, database_path = database_path,
                                       project_path = project_path){

    object <- format_chemical_properties(object)
    object <- tibble::tibble(object)

    readr::write_csv(x = object,
                     file = suppressWarnings(normalizePath(file.path(database_path, "chemical_properties.csv"))),
                     na = "NA")

    readr::write_csv(x = object,
                     file = suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))),
                     na = "NA")

}

# Save the project stage to RData object
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
