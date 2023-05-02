#' @importFrom purrr is_empty
#' @importFrom dplyr mutate left_join filter %>% select group_by ungroup
#' @importFrom dplyr pull arrange bind_rows rowwise if_else rename
#' @importFrom dplyr case_when slice_min summarise
#' @importFrom tibble tibble add_column add_row
#' @importFrom readr read_delim write_csv read_csv
#' @importFrom data.table %like% data.table
#' @importFrom progress progress_bar
#' @importFrom tidyr crossing
#' @importFrom webchem pc_prop
#' @importFrom utils globalVariables packageVersion
#' @importFrom EnvStats geoMean
#'
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\nThis is REcoTox version", packageVersion("REcoTox"), "\n"))
}


utils::globalVariables(
    names = c(
        ".", "..c2", ":=", "AVERAGE_MASS", "CASRN", "CID", "DTXSID_DTX", "EXCLUDE",
        "FOUND_BY", "INCHIKEY", "INCHI_STRING", "IUPAC_NAME", "LOG_S",
        "LOG_S_AD", "LOG_S_COMMENT", "MOLECULAR_FORMULA", "MONOISOTOPIC_MASS",
        "PREFERRED_NAME", "PubChem_CID", "QC_LEVEL", "QSAR_READY_SMILES",
        "REMARKS", "SMILES", "S_AD_index", "across", "author", "cas",
        "cas_number", "chemical_name", "cid", "common_name", "common_name_list",
        "compound_class", "conc1_comments", "conc1_max", "conc1_max_op",
        "conc1_mean", "conc1_mean_op", "conc1_min", "conc1_min_op", "conc1_type",
        "conc1_unit", "concentration_mean", "concentration_mean_list",
        "concentration_unit", "count_of_records", "database_path", "dtxsid",
        "dtxsid_ecotox", "duration_list", "ecotox_group", "effect",
        "effect_comments", "endpoint", "endpoint_comments",
        "endpoint_list", "family", "filename", "genus", "i.AVERAGE_MASS",
        "include_endpoint", "include_species", "input_data", "input_data_ad",
        "kingdom", "latin_name", "measurement", "measurement_comments",
        "measurement_list", "median", "obs_duration_mean", "obs_duration_unit",
        "organism_habitat", "organism_lifestage", "phylum_division", "project",
        "project_path", "publication_year", "reference_db", "reference_number",
        "reference_type", "result_id", "result_id_list", "species", "species_list",
        "species_number", "subphylum_div", "subspecies", "superclass",
        "tax_order", "test_cas", "test_id", "test_location", "title", "variety"
    ),
    package = "REcoTox", add = TRUE
)
