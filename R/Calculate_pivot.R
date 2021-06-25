calculate_pivot_table <- function(project, quantile = NA, limit_S_AD = NA, reread = FALSE){
  results <- data.table(project$object$results)
  if (reread == TRUE) {
    results <- fread(file = file.path(project_path,paste0(tolower(project$object$parameters$ecotoxgroup), "_final_results.csv")), na.strings = c("", NA), sep = ",", dec = ".", quote = "\"")
    
  }
  results_select <- results[which(results[, concentration_unit %like% "mg/L"])]
  results_select <- results_select[which(results_select[, concentration_unit > 0])]
  results_select <- results_select[!which(results_select[, EXCLUDE == 1])]
  #if(!is.na(limit_S_AD)){
  #  results_select <- results_select[which(results_select[,S_AD > limit_S_AD])]  
  #}
  #

  results_pivot <- results_select %>%  # Initial data
    tidyr::crossing(quantile = quantile) %>%  # Specify quantiles; crossing() is like expand.grid()
    dplyr::group_by(chemical_name, cas, cas_number, quantile) %>%  # Indicate your grouping var, plus your quantile var
    dplyr::summarise(quantile_value_mg_L = signif(quantile(concentration_mean, unique(quantile), na.rm =TRUE),4),
                     min_value_mg_L = signif(min(concentration_mean, na.rm = TRUE),4),
                     max_value_mg_L = signif(max(concentration_mean, na.rm = TRUE),4),
                     mean_value_mg_L = signif(mean(concentration_mean, na.rm = TRUE),4),
                     #sd_value_mg_L = signif(sd(concentration_mean, na.rm = TRUE),4),
                     geomean_value_mg_L = suppressWarnings(signif(EnvStats::geoMean(concentration_mean, na.rm = TRUE),4)),
                     #geosd_value_mg_L = suppressWarnings(signif(EnvStats::geoSD(concentration_mean, na.rm = TRUE),4)),
                     median_value_mg_L = signif(median(concentration_mean, na.rm = TRUE),4),
                     count = signif(sum(!is.na(concentration_mean)),4),
                     reference_count = length(unique(reference_number)),
                     #latin_name_count = length(unique(latin_name)),
                     OPERA_S_mg_L = signif(min(OPERA_S_mg_L, na.rm = TRUE), 4),
                     OPERA_AD_mean_S = mean(OPERA_S_mg_L_AD, na.rm = TRUE),
                     OPERA_AD_min_S = min(OPERA_S_mg_L_AD, na.rm = TRUE),
                     ACD_S_mg_L = signif(min(ACD_S_mg_L,na.rm = TRUE), 4),
                     ACD_AD_mean_S = min(ACD_S_mg_L_AD, na.rm = TRUE),
                     ACD_AD_min_S = min(ACD_S_mg_L_AD, na.rm = TRUE),
                     JC_S_mg_L = signif(min(JC_S_mg_L,na.rm = TRUE), 4),
                     JC_mean_AD_S = mean(JC_S_mg_L_AD, na.rm = TRUE),
                     JC_AD_min_S = min(JC_S_mg_L_AD, na.rm = TRUE),
                     Consensus_AD = signif(mean(c(OPERA_AD_mean_S,ACD_AD_mean_S, JC_mean_AD_S), na.rm = TRUE), 2),
                     duration_min_h = min(obs_duration_mean, na.rm = TRUE),
                     duration_max_h = max(obs_duration_mean, na.rm = TRUE),
                     
                     ) %>%  # unique() is needed
                    dplyr::mutate(quantile = sprintf("%1.0f%%", quantile * 100)) #%>% Optional prettification
                    

  # Add more tricky aggregated columns
  
  results_duration <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(duration_list = paste(obs_duration_mean, collapse = "|"))
  
  results_endpoints <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(endpoint_list = paste(endpoint, collapse = "|"))
  
  results_species <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(species_list = paste(latin_name, collapse = "|"))
  
  results_concentration <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(concentration_mean_list = paste(concentration_mean, collapse = "|"))
  
  results_measurement <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(measurement_list = paste(measurement, collapse = "|"))
  
  results_result_id <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(result_id_list = paste(result_id, collapse = "|"))
  
  results_common_name <- results_select %>%
    dplyr::group_by(chemical_name, cas, cas_number) %>%
    dplyr::summarise(common_name_list = paste(common_name, collapse = "|"))
  
  # Merge everything
  results_pivot <- left_join(results_pivot, results_duration[, c("duration_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_endpoints[, c("endpoint_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_species[, c("species_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_common_name[, c("common_name_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_result_id[, c("result_id_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_concentration[, c("concentration_mean_list", "cas_number")], by = "cas_number")
  results_pivot <- left_join(results_pivot, results_measurement[, c("measurement_list", "cas_number")], by = "cas_number")
  

  results_pivot <- left_join(results_pivot, project$object$chemprop[, c("cas_number", "DTXSID", "CASRN", "INCHIKEY", "PREFERRED_NAME",
                                                  "AVERAGE_MASS", "MONOISOTOPIC_MASS")], by = "cas_number")
  
  #Resample the output
  results_pivot <- results_pivot[c("cas_number", "PREFERRED_NAME", "quantile", "quantile_value_mg_L","min_value_mg_L",
                                   "max_value_mg_L",
                                   "mean_value_mg_L",
                                   #"sd_value_mg_L",
                                   "geomean_value_mg_L",
                                   #"geosd_value_mg_L",
                                   "median_value_mg_L",
                                   "count",
                                   "reference_count",
                                   #"latin_name_count",
                                   "Consensus_AD",
                                   "duration_min_h",
                                   "duration_max_h",
                                   "duration_list",
                                   "OPERA_S_mg_L",
                                   "ACD_S_mg_L",
                                   "JC_S_mg_L",
                                   "OPERA_AD_min_S",
                                   "ACD_AD_min_S",
                                   "JC_AD_min_S",
                                   "DTXSID",
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
                                   )]
  
  
  message("[EcoToxR]:  Save the pivot table.")
  fwrite(x = results_pivot, file = file.path(project_path,paste0(tolower(project$object$parameters$ecotoxgroup),"_EcoTox_pivot.csv")), na = NA, dec=".", sep=",")
  project$object$results_pivot <- results_pivot
  file <- suppressWarnings(normalizePath(file.path(project$project_path,paste0(tolower(project$object$parameters$ecotoxgroup),"_pivot_results.RData"))))
  message("[EcoToxR]:  Save the project.")
  save(project,file = file, compress = TRUE)
  return(project)
}
