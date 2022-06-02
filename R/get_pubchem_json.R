library(tidyverse)
library(jsonlite)
library(httr)


chemical_list <- read_csv("s:/Get_Pubchem/Basic_compound_list_20220109_original_DoNotChange_curated.csv", na = "n/a", show_col_types = FALSE)

chemical_list <- chemical_list %>% select(MoA_ID, PubChem_CID) %>% rename(cid = "PubChem_CID", ID = "MoA_ID")

#chemical_list_sample <- chemical_list %>% slice_sample(n = 25)

#chemical_list <- chemical_list[1:40,]


chemical_list_updated <- get_pubchem_experiments(object = chemical_list, sleep = 2, debug = FALSE)

get_pubchem_experiments <- function(object = chemical_list, sleep = 1, debug = FALSE){

    pubchem_experimental <- get_pubchem_experiments_tibble()

    length_progressbar <- nrow(object)
    pb <- progress::progress_bar$new(
        format = "[EcoToxR]:  Retrival of PubChem experimental data [:bar] :percent ETA: :eta",
        total = length_progressbar, clear = FALSE, width = 80)

    for ( i in 1:nrow(object) ){

        if (isTRUE(debug)) {
            print(i)
        }

        pb$tick()

        Sys.sleep(sleep)

        cid <- object[i, "cid"][[1]]
        ID <- object[i, "ID"][[1]]

        if ( is.na(cid) ) {

            pubchem_experimental_new_row <- get_pubchem_experiments_tibble() %>% add_row()

            pubchem_experimental_new_row$ID <- ID

            pubchem_experimental_new_row$cid <- cid

            pubchem_experimental <- rbind(pubchem_experimental, pubchem_experimental_new_row)

            next()

        } else {

            pubchem_json <- get_pubchem_json(cid = cid, ID = ID)


            if ( isTRUE(pubchem_json$json) ) {

                pubchem_json <- get_log_p(pubchem_json)

                pubchem_json <- get_ws(pubchem_json)

                pubchem_json <- get_pka(pubchem_json)

                pubchem_experimental_new_row <- get_pubchem_experiments_tibble() %>%
                    add_row() %>%
                    mutate(cid = pubchem_json$cid,
                           ID = pubchem_json$ID,
                           PC_LogP_experimental = pubchem_json$PC_LogP_experimental,
                           PC_LogP_experimental_reference = pubchem_json$PC_LogP_experimental_reference,
                           PC_WS_experimental = pubchem_json$PC_WS_experimental,
                           PC_WS_experimental_unit = pubchem_json$PC_WS_experimental_unit,
                           PC_WS_experimental_reference = pubchem_json$PC_WS_experimental_reference,
                           PC_pKa_experimental = pubchem_json$PC_pKa_experimental,
                           PC_pKa_experimental_reference = pubchem_json$PC_pKa_experimental_reference,
                           PC_acidic_pKa1_experimental = pubchem_json$PC_acidic_pKa1_experimental,
                           PC_acidic_pKa1_experimental_reference = pubchem_json$PC_acidic_pKa1_experimental_reference,
                           PC_acidic_pKa2_experimental = pubchem_json$PC_acidic_pKa2_experimental,
                           PC_acidic_pKa2_experimental_reference = pubchem_json$PC_acidic_pKa2_experimental_reference,
                           PC_acidic_pKa3_experimental = pubchem_json$PC_acidic_pKa3_experimental,
                           PC_acidic_pKa3_experimental_reference = pubchem_json$PC_acidic_pKa3_experimental_reference,
                           PC_basic_pKa1_experimental = pubchem_json$PC_basic_pKa1_experimental,
                           PC_basic_pKa1_experimental_reference = pubchem_json$PC_basic_pKa1_experimental_reference,
                           PC_basic_pKa2_experimental = pubchem_json$PC_basic_pKa2_experimental,
                           PC_basic_pKa2_experimental_reference = pubchem_json$PC_basic_pKa2_experimental_reference,
                           PC_basic_pKa3_experimental = pubchem_json$PC_basic_pKa3_experimental,
                           PC_basic_pKa3_experimental_reference = pubchem_json$PC_basic_pKa3_experimental_reference)

                pubchem_experimental <- rbind(pubchem_experimental, pubchem_experimental_new_row)
            }

            if ( isFALSE(pubchem_json$json) ) {

                pubchem_experimental_new_row <- get_pubchem_experiments_tibble() %>% add_row()

                pubchem_experimental_new_row$ID <- ID

                pubchem_experimental_new_row$cid <- cid

                pubchem_experimental <- rbind(pubchem_experimental, pubchem_experimental_new_row)

            }

        }

    }
    return(pubchem_experimental)
}



get_pubchem_json <- function(cid = NA, ID = NA) {

    pubchem_json <- list()
    pubchem_json$ID <- ID
    pubchem_json$cid <- cid

    url_json <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", cid, "/JSON")

    #resp <- tryCatch( {

    resp <- httr::GET(url_json)

    # },
    #
    # error = function(e) {
    #
    #   message("Server timed out")
    #
    #   resp <- 0
    #
    #   return(resp)
    #
    # }
    #
    # )
    #
    # if ( resp == 0 ) {
    #
    #   pubchem_json <- 0
    #
    #   break
    #
    #   return(pubchem_json)
    #
    # }
    #

    cont_raw <- tryCatch( {

        httr::content(resp, as = "parsed")

    },

    error = function(e) {

        message("JSON was empty")

        return()

    }
    )

    if ( is_empty(cont_raw) ) {

        pubchem_json$json <- FALSE

        return(pubchem_json)

    } else {

        pubchem_json$json <- TRUE

    }


    data_raw <- enframe(unlist(cont_raw))

    rgx_split <- "\\."
    n_cols_max <-
        data_raw %>%
        pull(name) %>%
        str_split(rgx_split) %>%
        map_dbl(~length(.)) %>%
        max()

    nms_sep <- paste0("name", 1:n_cols_max)
    data_sep <-
        data_raw %>%
        separate(name, into = nms_sep, sep = rgx_split, fill = "right")


    # Add column with row IDs

    data_sep <- data_sep %>% mutate(id = row_number())

    pubchem_json$data_sep <- data_sep

    pubchem_json$toc_heading_level3 <- data_sep %>%
        filter(
            name3 == "TOCHeading"
        )

    pubchem_json$toc_heading_level4 <- data_sep %>%
        filter(
            name4 == "TOCHeading"
        )

    pubchem_json$toc_heading_level5 <- data_sep %>%
        filter(
            name5 == "TOCHeading"
        )

    pubchem_json$unit_strings <- data_sep %>%
        filter(
            name7 == "Unit"
        )

    pubchem_json$string_markup <- data_sep %>%
        filter(
            name7 == "Markup"
        )

    pubchem_json$string_strings <- data_sep %>%
        filter(
            name7 == "String"
        )

    pubchem_json$numbers <- data_sep %>%
        filter(
            name7 == "Number"
        )


    pubchem_json$PC_LogP_experimental <- NA
    pubchem_json$PC_LogP_experimental_reference <- NA
    pubchem_json$PC_WS_experimental <- NA
    pubchem_json$PC_WS_experimental_reference <- NA
    pubchem_json$PC_WS_experimental_unit <- NA
    pubchem_json$PC_pKa_experimental <- NA
    pubchem_json$PC_pKa_experimental_reference <- NA
    pubchem_json$PC_pKa_experimental <- NA
    pubchem_json$PC_acidic_pKa1_experimental <- NA
    pubchem_json$PC_acidic_pKa1_experimental_reference <- NA
    pubchem_json$PC_acidic_pKa2_experimental <- NA
    pubchem_json$PC_acidic_pKa2_experimental_reference <- NA
    pubchem_json$PC_acidic_pKa3_experimental <- NA
    pubchem_json$PC_acidic_pKa3_experimental_reference <- NA
    pubchem_json$PC_basic_pKa1_experimental <- NA
    pubchem_json$PC_basic_pKa1_experimental_references <- NA
    pubchem_json$PC_basic_pKa2_experimental <- NA
    pubchem_json$PC_basic_pKa2_experimental_references <- NA
    pubchem_json$PC_basic_pKa3_experimental <- NA
    pubchem_json$PC_basic_pKa3_experimental_references <- NA


    return(pubchem_json)
}


get_log_p <- function(object = pubchem_json) {

    toc_level_5 <- object$toc_heading_level5
    data_sep <- object$data_sep

    if ( !is_empty(toc_level_5 %>% filter(value == "LogP") %>% pull(id)) ) {

        row_log_p <- toc_level_5 %>% filter(value == "LogP") %>% pull(id)

        row_log_p_row <- which(toc_level_5$id == row_log_p)

        row_log_p_1 <- toc_level_5 %>% slice(row_log_p_row + 1) %>% pull(id)

        log_p_section <- data_sep %>% slice(row_log_p:(row_log_p_1-1))

        number <- which(log_p_section$name7 == "Number")

        if ( is_empty(number) ) {

            string <- which(log_p_section$name8 == "String")
            object$PC_LogP_experimental <- log_p_section %>% slice(string[1]) %>% pull(value) %>% parse_number()
            object$PC_LogP_experimental_reference <- log_p_section %>% slice(string[1]-1) %>% pull(value)

        } else if ( is.na(number) ) {

            object$PC_LogP_experimental <- NA
            object$PC_LogP_experimental_reference <- NA

        } else {

            object$PC_LogP_experimental <- log_p_section %>% slice(number[1]) %>% pull(value) %>% parse_number()
            object$PC_LogP_experimental_reference <- log_p_section %>% slice(number[1]-1) %>% pull(value)

        }

    }

    return(object)
}

get_ws <- function(object = pubchem_json) {

    toc_level_5 <- object$toc_heading_level5
    data_sep <- object$data_sep

    if ( !is_empty(toc_level_5 %>% filter(value == "Solubility") %>% pull(id)) ) {

        row_s <- toc_level_5 %>% filter(value == "Solubility") %>% pull(id)
        row_s_row <- which(toc_level_5$id == row_s)

        row_s_1 <- toc_level_5 %>% slice(row_s_row + 1) %>% pull(id)

        s_section <- data_sep %>% slice(row_s:(row_s_1-1))

        number <- which(s_section$name7 == "Number")

        if ( !is_empty(number) ) {

            object$PC_WS_experimental <- s_section %>% slice(number[1]) %>% pull(value) %>% parse_number()

            object$PC_WS_experimental_reference <- s_section %>% slice(number[1] - 1) %>% pull(value)

            object$PC_WS_experimental_unit <- s_section %>% slice(number[1] + 1) %>% pull(value)

        }
    }

    return(object)
}


get_pka <- function(object = pubchem_json) {

    toc_level_5 <- object$toc_heading_level5

    data_sep <- object$data_sep

    if ( !is_empty(toc_level_5 %>% filter(value == "Dissociation Constants") %>% pull(id)) ) {

        row_dc <- toc_level_5 %>% filter(value == "Dissociation Constants") %>% pull(id)

        row_dc_row <- which(toc_level_5$id == row_dc)

        row_dc_1 <- pubchem_json$toc_heading_level5 %>% slice(row_dc_row + 1) %>% pull(id)

        dc_section <- pubchem_json$data_sep %>% slice(row_dc:(row_dc_1 - 1))

        pka <- which(dc_section$value == "pKa")

        acidic_pka <- which(dc_section$value == "Acidic pKa")

        basic_pka <- which(dc_section$value == "Basic pka")

        if ( is_empty(pka) & is_empty(acidic_pka) & is_empty(basic_pka) ) {

            pubchem_json$PC_pKa_experimental <- NA
            pubchem_json$PC_pKa_experimental_reference <- NA
            pubchem_json$PC_pKa_experimental <- NA
            pubchem_json$PC_acidic_pKa1_experimental <- NA
            pubchem_json$PC_acidic_pKa1_experimental_reference <- NA
            pubchem_json$PC_acidic_pKa2_experimental <- NA
            pubchem_json$PC_acidic_pKa2_experimental_reference <- NA
            pubchem_json$PC_acidic_pKa3_experimental <- NA
            pubchem_json$PC_acidic_pKa3_experimental_reference <- NA
            pubchem_json$PC_basic_pKa1_experimental <- NA
            pubchem_json$PC_basic_pKa1_experimental_references <- NA
            pubchem_json$PC_basic_pKa2_experimental <- NA
            pubchem_json$PC_basic_pKa2_experimental_references <- NA
            pubchem_json$PC_basic_pKa3_experimental <- NA
            pubchem_json$PC_basic_pKa3_experimental_references <- NA

        }


        if ( !is_empty(pka) ) {

            object$PC_pKa_experimental <- dc_section %>% slice(pka[1] + 2) %>% pull(value) %>% parse_number()

            object$PC_pKa_experimental_reference <- dc_section %>% slice(pka[1] + 1) %>% pull(value)

        }


        if ( !is_empty(acidic_pka) ) {

            object$PC_acidic_pKa1_experimental <- dc_section %>% slice(acidic_pka[1] + 3) %>% pull(value) %>% parse_number()

            object$PC_acidic_pKa1_experimental_reference <- dc_section %>% slice(acidic_pka[1] + 2) %>% pull(value)

            if ( length(acidic_pka) > 1 ) {

                object$PC_acidic_pKa2_experimental <- dc_section %>% slice(acidic_pka[2] + 3) %>% pull(value) %>% parse_number()

                object$PC_acidic_pKa2_experimental_reference <- dc_section %>% slice(acidic_pka[2] + 2) %>% pull(value)

            }

            if ( length(acidic_pka) > 2 ) {

                object$PC_acidic_pKa3_experimental <- dc_section %>% slice(acidic_pka[3] + 3) %>% pull(value) %>% parse_number()

                object$PC_acidic_pKa3_experimental_reference <- dc_section %>% slice(acidic_pka[3] + 2) %>% pull(value)

            }

        }


        if ( !is_empty(basic_pka) ) {

            object$PC_basic_pKa1_experimental <- dc_section %>% slice(basic_pka[1] + 3) %>% pull(value) %>% parse_number()

            object$PC_basic_pKa1_experimental_reference <- dc_section %>% slice(basic_pka[1] + 2) %>% pull(value)

            if ( length(basic_pka) > 1 ) {

                object$PC_basic_pKa2_experimental <- dc_section %>% slice(basic_pka[2] + 3) %>% pull(value) %>% parse_number()

                object$PC_basic_pKa1_experimental_reference <- dc_section %>% slice(basic_pka[2] + 2) %>% pull(value)

            }

            if ( length(basic_pka) > 2 ) {

                object$PC_basic_pKa3_experimental <- dc_section %>% slice(basic_pka[3] + 3) %>% pull(value) %>% parse_number()

                object$PC_basic_pKa3_experimental_reference <- dc_section %>% slice(basic_pka[3] + 2) %>% pull(value)

            }

        }
    }


    return(object)
}


get_pubchem_experiments_tibble <- function(){

    pubchem_experimental <- tibble(
        "ID" = character(),
        "cid" = integer(),
        "PC_LogP_experimental" = integer(),
        "PC_LogP_experimental_reference" = character(),
        "PC_WS_experimental" = integer(),
        "PC_WS_experimental_unit" = character(),
        "PC_WS_experimental_reference" = character(),
        "PC_pKa_experimental" = integer(),
        "PC_pKa_experimental_reference" = character(),
        "PC_acidic_pKa1_experimental" = integer(),
        "PC_acidic_pKa1_experimental_reference" = character(),
        "PC_acidic_pKa2_experimental" = integer(),
        "PC_acidic_pKa2_experimental_reference" = character(),
        "PC_acidic_pKa3_experimental" = integer(),
        "PC_acidic_pKa3_experimental_reference" = character(),
        "PC_basic_pKa1_experimental" = integer(),
        "PC_basic_pKa1_experimental_reference" = character(),
        "PC_basic_pKa2_experimental" = integer(),
        "PC_basic_pKa2_experimental_reference" = character(),
        "PC_basic_pKa3_experimental" = integer(),
        "PC_basic_pKa3_experimental_reference" = character()

    )

    return(pubchem_experimental)
}

