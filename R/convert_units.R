# Conversion functions

# Convert different volume units to mg/L or mol/L
convert_units <- function(object, sample_size = NA) {

    if(!is.na(sample_size)) {
        object <- object %>% dplyr::sample_n(size = sample_size)
    }

    message("[EcoToxR]:  Impute missing concentration mean by averaging minimum and maximum values")

    object <- object %>%
        dplyr::mutate(concentration_mean = conc1_mean, concentration_unit = conc1_unit)

    object <-
        dplyr::bind_rows(object %>%
                             dplyr::filter(!is.na(concentration_mean)),
                         object %>%
                             dplyr::filter(is.na(concentration_mean)) %>%
                             dplyr::rowwise() %>%
                             dplyr::mutate(concentration_mean =
                                               dplyr::if_else(condition = is.na(concentration_mean),
                                                              true = mean(c(conc1_min, conc1_max)),
                                                              false = concentration_mean))
        )


    # convert similar unit to SI conform units (e.g. ppm to mg/L)

    # ng related
    message("[EcoToxR]:  Converting units to human readible format (i.e. mg/L).")
    message("[EcoToxR]:  The following units will be converted:")
    print(unique(object$concentration_unit))

    # this assumes that only some records are reported in log(LC50)
    log <- "^.log.[A-Z]*[0-9]*?$"

    if (nrow(object %>%
             dplyr::filter(endpoint %like% log)) > 0) {

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = endpoint %like% log,
                                             true = 10^concentration_mean,
                                             false= concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(endpoint =
                              dplyr::if_else(condition = concentration_unit %like% log,
                                             true = "LC50",
                                             false = endpoint))
    }

    # pg/L
    pg <- "^pg/(l|L)$"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% pg)) > 0) {

        message("[EcoToxR]:  Converting pg/L like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% pg,
                                             true = concentration_mean / 1e+09,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% pg,
                                             true = "mg/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping pg/L like units")}


    #ng
    ng <- "^(|A(I|E|i|e) )ng/(l|L|dm3)$|^ppt$|^(|A(I|E|i|e) )pg/(ml|mL)$"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% ng)) > 0) {

        message("[EcoToxR]:  Converting ng/L like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% ng,
                                             true = concentration_mean / 1e+06,
                                             false = concentration_mean))
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% ng,
                                             true ="mg/L",
                                             false =concentration_unit))

    } else {message("[EcoToxR]:  Skipping ng/L like units")}

    #µg related
    ug <- "(^(|(A|a)(I|E|i|e) )ug/(l|L|dm3)$)|^mg/m3$|^ppb$|(^pg/u(l|L)$)|(^(|(A|a)(I|E|i|e) )ng/m(l|L)$)|^pg/u(l|L)"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% ug)) > 0) {

        message("[EcoToxR]:  Converting ug/L like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% ug,
                                             true = concentration_mean / 1e+03,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% ug,
                                             true = "mg/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping ug/L like units")}


    # g/L
    g <- "(^(|(A|a)(e|i|E|I) )(g/(l|L)$))|(^mg/(ml|mL)$|(^g/dm3$))|(^ug/(ul|uL)$)"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% g)) > 0) {

        message("[EcoToxR]:  Converting g/L like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% g,
                                             true = concentration_mean * 1e+03,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% g,
                                             true = "mg/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping g/L like units")}


    # mg/L
    message("[EcoToxR]:  Converting mg/L like units")
    mg <- "(^(|((A|a)(I|E|i|e)) )mg/(l|L|dm3)$)|^ppm$|(^(|((A|a)(I|E|i|e)) )ug/m(l|L|m3)$)|^g/m3$"
    object <-
        object %>%
        dplyr::rowwise() %>%
        dplyr::mutate(concentration_unit =
                          dplyr::if_else(condition = concentration_unit %like% mg,
                                         true = "mg/L",
                                         false = concentration_unit))

    #pmol related
    pmol <- "(^(|(A|a)(I|E|i|e) )pmol/(l|L)$)"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% pmol)) > 0) {

        message("[EcoToxR]:  Converting pmol like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% pmol,
                                             true = concentration_mean / 1e+12,
                                             false = concentration_mean))
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% pmol,
                                             true = "mol/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping pmol/L like units")}

    #nmol related
    nmol <- "(^(|(A|a)(I|E|i|e) )nmol/(l|L)$)"

    if (nrow(object %>%
             dplyr::filter(concentration_unit %like% nmol)) > 0) {

        message("[EcoToxR]:  Converting nmol like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% nmol,
                                             true = concentration_mean / 1e+09,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% nmol,
                                             true = "mol/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping nmol/L like units")}

    #umol related
    umol <- "(^(|(A|a)(I|E|i|e) )umol/(dm3|l|L)$)|^nmol/ml$|^mmol/m3$"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% umol)) > 0) {

        message("[EcoToxR]:  Converting pmol like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% umol,
                                             true = concentration_mean / 1e+06,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% umol,
                                             true = "mol/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping umol/L like units")}

    #mmol related
    mmol <- "^(|(A|a)(I|E|i|e) )mmol/(dm3|l|L)$"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% mmol)) > 0) {

        message("[EcoToxR]:  Converting mmol like units")
        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_mean =
                              dplyr::if_else(condition = concentration_unit %like% mmol,
                                             true = concentration_mean / 1e+03,
                                             false = concentration_mean))

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% mmol,
                                             true = "mol/L",
                                             false = concentration_unit))

    } else {message("[EcoToxR]:  Skipping mmol/L like units")}


    #mol related
    mol <- "^(|(A|a)(I|E|i|e) )mol/(dm3)$"

    if (nrow(object %>% dplyr::filter(concentration_unit %like% mol)) > 0) {
        message("[EcoToxR]:  Converting mol like units")

        object <-
            object %>%
            dplyr::rowwise() %>%
            dplyr::mutate(concentration_unit =
                              dplyr::if_else(condition = concentration_unit %like% mol,
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

# Convert the water solubility in log mol/L to mg/L
convert_water_solubility <- function(object) {

    message("[EcoToxR]:  Convert the solubility to mg/L.")

    object <- object %>%
        dplyr::rowwise() %>%
        dplyr::mutate(S_mg_L = dplyr::if_else(!is.na(AVERAGE_MASS),
                                              true = signif(x = 10^LOG_S * 1000 * AVERAGE_MASS, digits = 4),
                                              false = NaN)) %>%
        ungroup() %>%
        dplyr::rename(QSAR_S_AD = "LOG_S_AD", QSAR_S_COMMENT = "LOG_S_COMMENT")

    return(object)

}

# Update the remaining values in mol/L to mg/L after review of the mol table
update_mol_units <- function(object, database_path = project$database_path, project_path = project$project_path) {
    message("[EcoToxR]:  Recalulating record unit mol/L to mg/L.")
    ecotoxgroup <- object$parameters$ecotoxgroup
    file_name <- file.path(project_path, paste0(tolower(ecotoxgroup), "_mol_weight.csv"))

    chemprop <- readr::read_csv(file = file_name,
                                col_names = TRUE,
                                na = c("NA", "", NA, NaN, "N/A", "n/a"),
                                show_col_types = FALSE)

    chemprop <- format_chemical_properties(chemprop)

    # subset

    mol_records <- object$results_filtered %>%
        dplyr::filter(concentration_unit %like% "mol/L") %>%
        dplyr::left_join(chemprop %>%
                             dplyr::select(cas_number, AVERAGE_MASS),
                         by = "cas_number") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(concentration_mean = concentration_mean * 1000 * AVERAGE_MASS) %>%
        dplyr::mutate(concentration_unit = "mg/L") %>%
        dplyr::select(-AVERAGE_MASS)

    object$results_filtered <-
        dplyr::bind_rows(object$results_filtered %>%
                             dplyr::filter(concentration_unit == "mg/L"), mol_records)

    message("[EcoToxR]: The data finally contains the following units ")
    print(unique(object$results_filtered$concentration_unit))


    object$chemprop <-
        tibble::tibble(data.table(object$chemprop)[data.table(chemprop),
                                                   on = c("cas_number"),
                                                   AVERAGE_MASS := i.AVERAGE_MASS])

    readr::write_csv(x = object$chemprop,
                     file = suppressWarnings(normalizePath(file.path(project_path, "chemical_properties.csv"))),
                     col_names = TRUE)

    return(object)
}
