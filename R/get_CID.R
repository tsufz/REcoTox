library(webchem)
library(progress)
library(tidyverse)
library(data.table)

chemical_list <- read_csv("c:/TEMP/EcoToxDB/test/algae_chemical_list.csv")

length_progressbar <- nrow(chemical_list)
pb <- progress::progress_bar$new(
    format = "[EcoToxR]:  Retrival of PubChem data [:bar] :percent ETA: :eta",
    total = length_progressbar, clear = FALSE, width = 80)

for (i in 1:nrow(chemical_list)){

    pb$tick()

    if (is.na(chemical_list$CID[i])){
    # lookup for CID based on CASRN
    casrn <- chemical_list[i, "CASRN"][[1]]
    pccid <- get_cid(casrn)
    cas_number <- chemical_list[i, "cas_number"][[1]]
    pccid <- pccid %>% mutate(across(cid, as.integer)) %>% mutate(cas_number = cas_number)


    if (is.na(pccid$cid[[1]])) {
        next()

    } else if (nrow(pccid) > 1) {
        # Us the first entry in PubChem only
        pccid <- pccid %>% arrange(cid)
        pccid <- pccid %>% slice_min(cid, n = 1)

    }


    chemical_list[i, "CID"] <- pccid[[2]]

    # write_csv(chemical_list, "c:/TEMP/EcoToxDB/test/algae_chemical_list.csv")

    } else {
        next()
    }

}



write_csv(chemical_list, "c:/TEMP/EcoToxDB/test/algae_chemical_list_update.csv")
