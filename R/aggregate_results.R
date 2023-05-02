# Workflow to aggregate the longer tabulated results of the
# EcoTox Knowledgebase query to wider format. The table is exported
# in \code{csv} format to the \code{project_folder}.
#
#' @title Aggregate data
#'
#' @description
#' This function aggregates the longer formmat results of the EcoTox
#' Knowledgebase query stored in the \code{REcoTox} project to the
#' final longer results data table.
#'
#' @param project Name of the initial project in the \code{environment}.
#' The default value is \code{project}.
#'
#' @param quantile Value of the lowest quantile to aggregate the data.
#' Default value: \code{NA}. For example: \code{0.05}
#'
#' @param file_name Custom file name (without extension) for the exported
#' results. Default value: \code{NA}.
#'
#' @param reread Read the longer \code{csv} table from \code{project_folder}
#' and process it. For example to exclude compounds or single results.
#' Default value: \code{FALSE} (values: \code{c(FALSE, TRUE)}).
#'
#' @param save_project_steps Stores the current project state in the \code{project_folder}.
#' The default value is \code{FALSE} (values: \code{c(FALSE, TRUE)}).
#'
#' @author Tobias Schulze
#'
#' @examples
#' # Run the current project with a quantile of {0.05} with a custom name.
#'
#' \dontrun{aggregate_results(project = project, quantile = 0.05,
#' file_name = "algae_aggregated_results")}
#'
#' @export
#'

aggregate_results <- function(project = project,
                                  quantile = NA,
                                  file_name = NA,
                                  reread = FALSE,
                                  save_project_steps = FALSE) {

message("[EcoToxR]:  Summarizing the data.")

  results <- project$object$results
  chemical_list <- project$object$chemical_list

  if (reread == TRUE) {
    results <- readr::read_csv(file = file.path(project_path, paste0(tolower(project$object$parameters$ecotoxgroup), "_final_results.csv")), na = c("", NA), show_col_types = FALSE)

  }

  results <- results %>%
    dplyr::group_by(cas_number) %>%
    dplyr::arrange(.by_group = TRUE)

  results_select <- results %>%
    dplyr::filter(concentration_unit == "mg/L", concentration_mean > 0, is.na(EXCLUDE))

  results_pivot <- results_select %>%  # Initial data
    tidyr::crossing(quantile = quantile) %>%  # Specify quantiles; crossing() is like expand.grid()
    dplyr::group_by(cas_number) %>%  # Indicate your grouping var, plus your quantile var
    dplyr::summarise(quantile_value_mg_L = signif(quantile(concentration_mean, unique(quantile), na.rm = TRUE), 4),
              min_value_mg_L = signif(min(concentration_mean, na.rm = TRUE), 4),
              max_value_mg_L = signif(max(concentration_mean, na.rm = TRUE), 4),
              mean_value_mg_L = signif(mean(concentration_mean, na.rm = TRUE), 4),
              #sd_value_mg_L = signif(sd(concentration_mean, na.rm = TRUE),4),
              geomean_value_mg_L = suppressWarnings(signif(EnvStats::geoMean(concentration_mean, na.rm = TRUE), 4)),
              #geosd_value_mg_L = suppressWarnings(signif(EnvStats::geoSD(concentration_mean, na.rm = TRUE),4)),
              median_value_mg_L = signif(median(concentration_mean, na.rm = TRUE), 4),
              result_count = signif(sum(!is.na(concentration_mean)), 4),
              reference_count = length(unique(reference_number)),
              #latin_name_count = length(unique(latin_name)),
              duration_min_h = min(obs_duration_mean, na.rm = TRUE),
              duration_max_h = max(obs_duration_mean, na.rm = TRUE),
              .groups = "keep"
    ) %>%
    dplyr::ungroup()

  # Add more tricky aggregated columns

  results_duration <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(duration_list = paste(obs_duration_mean, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, duration_list)

  results_endpoints <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(endpoint_list = paste(endpoint, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, endpoint_list)

  results_species <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr:: summarise(species_list = paste(latin_name, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, species_list)

  results_concentration <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(concentration_mean_list = paste(concentration_mean, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, concentration_mean_list)

  results_measurement <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(measurement_list = paste(measurement, collapse = "|"), .groups = "keep")  %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, measurement_list)

  results_result_id <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(result_id_list = paste(result_id, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, result_id_list)

  results_common_name <- results_select %>%
    dplyr::group_by(cas_number) %>%
    dplyr::summarise(common_name_list = paste(common_name, collapse = "|"), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::select(cas_number, common_name_list)

  # Merge everything
  results_pivot <- results_pivot %>%
    dplyr::left_join(results_duration %>% dplyr::select(duration_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_endpoints %>% dplyr::select(endpoint_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_species %>% dplyr::select(species_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_common_name %>% dplyr::select(common_name_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_result_id %>% dplyr::select(result_id_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_concentration %>% dplyr::select(concentration_mean_list, cas_number), by = "cas_number") %>%
    dplyr::left_join(results_measurement %>% dplyr::select(measurement_list, cas_number), by = "cas_number")

  results_pivot <-
    results_pivot %>%
    dplyr::left_join(chemical_list %>%
                       dplyr::select(cas_number,
                                     chemical_name,
                                     cas,
                                     dtxsid_ecotox,
                                     LOG_S,
                                     LOG_S_AD,
                                     S_AD_index,
                                     LOG_S_COMMENT,
                                     AVERAGE_MASS,
                                     DTXSID_DTX,
                                     PubChem_CID,
                                     CASRN,
                                     SMILES,
                                     INCHIKEY,
                                     PREFERRED_NAME,
                                     IUPAC_NAME,
                                     INCHIKEY,
                                     QSAR_READY_SMILES,
                                     MOLECULAR_FORMULA,
                                     REMARKS),
                                     by = "cas_number") %>%
    dplyr::mutate("quantile" = quantile)

  results_pivot <- convert_water_solubility(object = results_pivot)

  results_pivot <- calculate_solubility_domain(object = results_pivot)

  #Resample the output
  results_pivot <- results_pivot %>% dplyr::select("cas_number",
                                   "chemical_name",
                                   "PREFERRED_NAME",
                                   "DTXSID_DTX",
                                   "PubChem_CID",
                                   "quantile",
                                   "quantile_value_mg_L",
                                   "min_value_mg_L",
                                   "max_value_mg_L",
                                   "mean_value_mg_L",
                                   "geomean_value_mg_L",
                                   "median_value_mg_L",
                                   "quantile_value_S_AD",
                                   "min_value_S_AD",
                                   "max_value_S_AD",
                                   "mean_value_S_AD",
                                   "geomean_value_S_AD",
                                   "median_value_S_AD",
                                   "S_mg_L",
                                   "QSAR_S_AD",
                                   "S_AD_index",
                                   "QSAR_S_COMMENT",
                                   "duration_min_h",
                                   "duration_max_h",
                                   "result_count",
                                   "reference_count",
                                   "IUPAC_NAME",
                                   "CASRN",
                                   "SMILES",
                                   "INCHIKEY",
                                   "QSAR_READY_SMILES",
                                   "MOLECULAR_FORMULA",
                                   "AVERAGE_MASS",
                                   "REMARKS",
                                   "endpoint_list",
                                   "species_list",
                                   "concentration_mean_list",
                                   "measurement_list",
                                   "result_id_list",
                                   "duration_list"
                                   ) %>%
    dplyr::group_by(chemical_name) %>%
    dplyr::arrange(.by_group = TRUE) %>%
    dplyr::ungroup()


  message("[EcoToxR]:  Saving the aggregated table.")

  if (is.na(file_name)) {
    file = file.path(project_path,
                     paste0(tolower(project$object$parameters$ecotoxgroup),
                            "_EcoToxKB_aggregated_results.csv"))
  } else {
    file = file.path(project_path,
                     paste0(tolower(filename),
                            ".csv"))
  }

  readr::write_csv(x = results_pivot, file = file, na = "NA", col_names = TRUE)
  project$object$results_pivot <- results_pivot

  if (isTRUE(save_project_steps)) {
    file <- suppressWarnings(normalizePath(file.path(project$project_path, paste0(tolower(project$object$parameters$ecotoxgroup), "_aggregated_results.RData"))))
    message("[EcoToxR]:  Saving the project.")
    save(project,file = file, compress = TRUE)
  }
  return(project)
}
