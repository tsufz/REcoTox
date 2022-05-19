calculate_pivot_table <- function(project, quantile = NA, limit_S_AD = NA, reread = FALSE){
  results <- project$object$results
  if (reread == TRUE) {
    results <- read_csv(file = file.path(project_path, paste0(tolower(project$object$parameters$ecotoxgroup), "_final_results.csv")), na = c("", NA), show_col_types = FALSE)

  }

  results_select <- results %>% filter(concentration_unit == "mg/L", concentration_mean > 0, is.na(EXCLUDE))


  #results_select <- results[which(results[, concentration_unit %like% "mg/L"])]
  #results_select <- results_select[which(results_select[, concentration_unit > 0])]
  #results_select <- results_select[!which(results_select[, EXCLUDE == 1])]

  results_pivot <- results_select %>%  # Initial data
    crossing(quantile = quantile) %>%  # Specify quantiles; crossing() is like expand.grid()
    group_by(chemical_name, cas, cas_number, quantile) %>%  # Indicate your grouping var, plus your quantile var
          summarise(quantile_value_mg_L = signif(quantile(concentration_mean, unique(quantile), na.rm = TRUE), 4),
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
                    S_mg_L =  suppressWarnings(if_else(condition = !is.infinite(min(S_mg_L, na.rm = TRUE)),
                                                           true = signif(min(S_mg_L, na.rm = TRUE), 4),
                                                           false = NaN)),
                    AD_mean_S =  suppressWarnings(if_else(condition = !is.infinite(mean(S_mg_L_AD, na.rm = TRUE)),
                                                              true = signif(mean(S_mg_L_AD, na.rm = TRUE), 1),
                                                              false = NaN)),
                    AD_min_S =  suppressWarnings(if_else(condition = !is.infinite(min(S_mg_L_AD, na.rm = TRUE)),
                                                             true = signif(min(S_mg_L_AD, na.rm = TRUE), 1),
                                                             false = NaN)),
                    duration_min_h = min(obs_duration_mean, na.rm = TRUE),
                    duration_max_h = max(obs_duration_mean, na.rm = TRUE)
          )

                    #%>%  # unique() is needed
                    #   dplyr::mutate(quantile = sprintf("%1.0f%%", quantile * 100)) #%>% Optional prettification


  # Add more tricky aggregated columns

  results_duration <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(duration_list = paste(obs_duration_mean, collapse = "|"))

  results_endpoints <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(endpoint_list = paste(endpoint, collapse = "|"))

  results_species <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(species_list = paste(latin_name, collapse = "|"))

  results_concentration <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(concentration_mean_list = paste(concentration_mean, collapse = "|"))

  results_measurement <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(measurement_list = paste(measurement, collapse = "|"))

  results_result_id <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(result_id_list = paste(result_id, collapse = "|"))

  results_common_name <- results_select %>%
    group_by(chemical_name, cas, cas_number) %>%
    summarise(common_name_list = paste(common_name, collapse = "|"))

  # Merge everything
  results_pivot <- results_pivot %>%
      left_join(results_duration %>% select(duration_list, cas_number), by = "cas_number") %>%
      left_join(results_endpoints %>% select(endpoint_list, cas_number), by = "cas_number") %>%
      left_join(results_species %>% select(species_list, cas_number), by = "cas_number") %>%
      left_join(results_common_name %>% select(common_name_list, cas_number), by = "cas_number") %>%
      left_join(results_result_id %>% select(result_id_list, cas_number), by = "cas_number") %>%
      left_join(results_concentration %>% select(concentration_mean_list, cas_number), by = "cas_number") %>%
      left_join(results_measurement %>% select(measurement_list, cas_number), by = "cas_number")


  results_pivot <- results_pivot %>%
      left_join(project$object$chemprop %>% select("cas_number", "CID", "DTXSID", "CASRN", "INCHIKEY", "PREFERRED_NAME",
                                                  "AVERAGE_MASS", "MONOISOTOPIC_MASS"), by = "cas_number")

  #Resample the output
  results_pivot <- results_pivot %>% select("cas_number",
                                   "PREFERRED_NAME",
                                   "quantile",
                                   "quantile_value_mg_L",
                                   "min_value_mg_L",
                                   "max_value_mg_L",
                                   "mean_value_mg_L",
                                   #"sd_value_mg_L",
                                   "geomean_value_mg_L",
                                   #"geosd_value_mg_L",
                                   "median_value_mg_L",
                                   "result_count",
                                   "reference_count",
                                   #"latin_name_count",
                                   "duration_min_h",
                                   "duration_max_h",
                                   "duration_list",
                                   "S_mg_L",
                                   "AD_mean_S",
                                   "DTXSID",
                                   "CID",
                                   "CASRN",
                                   "INCHIKEY",
                                   "AVERAGE_MASS",
                                   "endpoint_list",
                                   "species_list",
                                   "concentration_mean_list",
                                   "measurement_list",
                                   "result_id_list"



                                   #"value_list",
                                   #"endpoint_list",
                                   #"measurement_list",
                                   #"species_list",
                                   #"test_list"
                                   ) %>%
      group_by(PREFERRED_NAME) %>%
      arrange(.)


  message("[EcoToxR]:  Saving the pivot table.")
  write_csv(x = results_pivot, file = file.path(project_path,paste0(tolower(project$object$parameters$ecotoxgroup), "_EcoTox_pivot.csv")), na = "NA", col_names = TRUE)
  project$object$results_pivot <- results_pivot
  file <- suppressWarnings(normalizePath(file.path(project$project_path,paste0(tolower(project$object$parameters$ecotoxgroup), "_pivot_results.RData"))))
  message("[EcoToxR]:  Saving the project.")
  save(project,file = file, compress = TRUE)
  return(project)
}
