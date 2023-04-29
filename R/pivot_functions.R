calculate_hours <- function(object = object) {
    message("[EcoToxR]:  Recalculating duration data (day to hours).")
    results <- tibble::tibble(object$results)

    results_h <- results %>% dplyr::filter(obs_duration_unit == "h")
    results_d <- results %>% dplyr::filter(obs_duration_unit == "d")

    results_d <- results_d %>%
        dplyr::filter(obs_duration_unit == "d") %>%
        dplyr::mutate(obs_duration_mean = 24 * obs_duration_mean) %>%
        dplyr::mutate(obs_duration_unit = "h")

    results <- tibble::tibble(dplyr::bind_rows(results_h, results_d))

    object$results <- results

    return(object)

}

calculate_solubility_domain <-
    function(object = object, input_column_list = c("quantile_value_mg_L",
                                                    "min_value_mg_L",
                                                    "max_value_mg_L",
                                                    "mean_value_mg_L",
                                                    "geomean_value_mg_L",
                                                    "median_value_mg_L")) {

        message("[EcoToxR]:  Estimating the solubility domain.")

        for (input_column in input_column_list) {

            object <- object %>% ungroup()

            object <- object %>%
                dplyr::mutate(input_data =
                                  as.vector(dplyr::pull(object[, input_column]))) %>%
                tibble::add_column(input_data_ad = NA) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(input_data_ad =
                                  dplyr::case_when(

                                      input_data <= S_mg_L ~ 3,

                                      input_data > S_mg_L & input_data <= 10^5 * log10(S_mg_L) ~ 2,

                                      input_data > S_mg_L & input_data > 10^5 * log10(S_mg_L) &
                                          input_data <= 10^10 * log10(S_mg_L) ~ 1,

                                      input_data > 10^10 * log10(S_mg_L) ~ 0
                                  )
                )


            ad_output_column =
                dplyr::case_when(input_column == "quantile_value_mg_L" ~ "quantile_value_S_AD",
                                 input_column == "min_value_mg_L" ~ "min_value_S_AD",
                                 input_column == "max_value_mg_L" ~ "max_value_S_AD",
                                 input_column == "mean_value_mg_L" ~ "mean_value_S_AD",
                                 input_column == "geomean_value_mg_L" ~ "geomean_value_S_AD",
                                 input_column == "median_value_mg_L" ~ "median_value_S_AD")


            object <- object %>%
                dplyr::select(., -input_data) %>%
                dplyr::rename(!!ad_output_column := input_data_ad) %>%
                dplyr::ungroup()

        }

        return(object)
    }
build_final_list <- function(object = object){

    results <- object$results_filtered

    if(is.null(results$include_species)){
        results[, "include_species" := NA]
    }


    results <- results %>% dplyr::select(cas_number, cas, chemical_name,
                                         compound_class, dtxsid_ecotox,
                                         test_location, conc1_type, conc1_mean_op,
                                         conc1_mean, conc1_min_op, conc1_min,
                                         conc1_max_op, conc1_max, conc1_unit,
                                         conc1_comments, concentration_mean,
                                         concentration_unit, obs_duration_mean,
                                         obs_duration_unit, test_id,
                                         reference_number, endpoint,
                                         endpoint_comments, effect,
                                         effect_comments, measurement,
                                         measurement_comments, organism_lifestage,
                                         common_name, latin_name, kingdom,
                                         phylum_division, subphylum_div,
                                         superclass, class, tax_order, family,
                                         genus, species, subspecies,
                                         variety, ecotox_group, result_id,
                                         reference_db, reference_type, author,
                                         title, source, publication_year,
                                         include_endpoint, include_species) %>%
        dplyr::rename(dtxsid_ecotox = dtxsid_ecotox)


    results <- results %>%
        dplyr::left_join(object$chemprop %>%
                             dplyr::select(cas_number,
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
                                           S_AD_index,
                                           LOG_S_COMMENT,
                                           EXCLUDE,
                                           REMARKS),
                         by = "cas_number"
        )


    # reorder list

    results <- results %>% dplyr::select(EXCLUDE, REMARKS, cas_number, cas,
                                         PubChem_CID, dtxsid_ecotox, DTXSID_DTX,
                                         PREFERRED_NAME, test_location,
                                         reference_number, conc1_type,
                                         conc1_mean_op, conc1_mean, conc1_min_op,
                                         conc1_min, conc1_max_op, conc1_max,
                                         conc1_unit, conc1_comments,
                                         concentration_mean, concentration_unit,
                                         test_id, endpoint, endpoint_comments,
                                         effect, effect_comments, measurement,
                                         measurement_comments, obs_duration_mean,
                                         obs_duration_unit, organism_lifestage,
                                         common_name, latin_name, kingdom,
                                         phylum_division, subphylum_div,
                                         superclass, class, tax_order, family,
                                         genus, species, subspecies, variety,
                                         ecotox_group, chemical_name,
                                         compound_class, CASRN, SMILES, INCHIKEY,
                                         QSAR_READY_SMILES, MOLECULAR_FORMULA,
                                         AVERAGE_MASS, MONOISOTOPIC_MASS, LOG_S,
                                         LOG_S_AD, S_AD_index, LOG_S_COMMENT,
                                         reference_db, result_id,
                                         reference_type, author, title, source,
                                         publication_year, include_endpoint,
                                         include_species)
    object$results <- results
    return(object)
}
