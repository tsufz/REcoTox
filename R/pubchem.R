library(webchem)
library(progress)
library(tidyverse)
library(data.table)

chemical_list <- read_csv("C:/TEMP/EcoToxDB/test/pubchem.csv", na = c("N/A", "NA", NA))

# Split the list in entries with smiles and w/o smiles
chemical_list_SMILES <- chemical_list %>% filter(!is.na(SMILES))
chemical_list_no_SMILES <- chemical_list %>% filter(is.na(SMILES))




# get information from pubchem to fill gaps of DTXSID query
#
#
#
pubchem <- tibble(
    "cas_number" = integer(),
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
    pccid <- get_cid(casrn)
    cas_number <- chemical_list_no_SMILES[i, "cas_number"][[1]]
    pccid <- pccid %>% mutate(across(cid, as.integer)) %>% mutate(cas_number = cas_number)


    if (is.na(pccid$cid[[1]])) {
        pubchem_new_row <- tibble(
            "cas_number" = integer(),
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

        pubchem_new_row <- pubchem_new_row %>% add_row() %>% mutate(cas_number = cas_numb, cas = casrn, FOUND_BY = "No data retrieved from PubChem")

        pubchem <- rbind(pubchem, pubchem_new_row)

        next()

    } else if (nrow(pccid) > 1) {
        # Us the first entry in PubChem only
        pccid <- pccid %>% arrange(cid)
        pccid <- pccid %>% slice_min(cid, n = 1)

    } else {
        pc_props <- tibble(pc_prop(pccid$cid, properties = c("Title",
                                                              "InChIKey",
                                                              "IUPACName",
                                                              "CanonicalSMILES",
                                                              "InChI",
                                                              "MolecularFormula",
                                                              "MolecularWeight",
                                                              "MonoisotopicMass")))

        if (is.na(pc_props$CanonicalSMILES)) {

            pubchem_new_row <- tibble(
                "cas_number" = integer(),
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

            pubchem_new_row <- pubchem_new_row %>% add_row() %>% mutate(cas_number = cas_numb, cas = casrn, FOUND_BY = "No data retrieved from PubChem")

            pubchem <- rbind(pubchem, pubchem_new_row)

            next()
        }


        # Postprocess the retrieved data
        pc_props <- pccid %>% left_join(pc_props, by = c("cid" = "CID"))

        if (!is_empty(which(names(pc_props) %like% "IUPACName"))) {

            pc_props <- pc_props %>% rename(c("cas" = "query", "cid" = "cid", "cas_number" = "cas_number",
                                          "PREFERRED_NAME" = "Title",
                                          "MOLECULAR_FORMULA" = "MolecularFormula",
                                          "AVERAGE_MASS" = "MolecularWeight", "SMILES" = "CanonicalSMILES",
                                          "INCHI_STRING"  = "InChI", "INCHIKEY" = "InChIKey",
                                          "IUPAC_NAME" = "IUPACName", "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

            pc_props <- pc_props %>% mutate("CASRN" = cas)

            pc_props <- pc_props %>% select(cas_number, cas, PREFERRED_NAME, CASRN,
                                            INCHIKEY, IUPAC_NAME, SMILES, INCHI_STRING, MOLECULAR_FORMULA,
                                            AVERAGE_MASS, MONOISOTOPIC_MASS)

        } else {

            pc_props <- pc_props %>% rename(c("cas" = "query", "cid" = "cid", "cas_number" = "cas_number",
                                            "PREFERRED_NAME" = "Title",
                                            "MOLECULAR_FORMULA" = "MolecularFormula",
                                            "AVERAGE_MASS" = "MolecularWeight", "SMILES" = "CanonicalSMILES",
                                            "INCHI_STRING"  = "InChI", "INCHIKEY" = "InChIKey",
                                            "MONOISOTOPIC_MASS" = "MonoisotopicMass"))

            pc_props <- pc_props %>% mutate(IUPAC_NAME = NA)
            pc_props <- pc_props %>% mutate(CASRN = cas)

            pc_props <- pc_props %>% select(cas_number, cas, PREFERRED_NAME, CASRN,
                                            INCHIKEY, IUPAC_NAME, SMILES, INCHI_STRING, MOLECULAR_FORMULA,
                                            AVERAGE_MASS, MONOISOTOPIC_MASS)
        }

        pc_props <- pc_props %>% add_column(FOUND_BY = "PubChem", .before = 3)


    }

    pubchem <- pubchem %>% rbind(pc_props)

}

# Output of webchem is character, needs to be fixed here.
pubchem <- pubchem %>% mutate(across(cas_number, as.integer)) %>% mutate(across(AVERAGE_MASS:MONOISOTOPIC_MASS, as.double))

chemicals_update <- chemical_list_no_SMILES %>% inner_join(pubchem %>% select(cas_number))

chemicals_update <- chemicals_update %>%
    mutate(FOUND_BY = pubchem$FOUND_BY, PREFERRED_NAME = pubchem$PREFERRED_NAME, CASRN = pubchem$CASRN,
           INCHIKEY = pubchem$INCHIKEY, IUPAC_NAME = pubchem$IUPAC_NAME, SMILES = pubchem$SMILES,
           INCHI_STRING = pubchem$INCHI_STRING, MOLECULAR_FORMULA = pubchem$MOLECULAR_FORMULA,
           AVERAGE_MASS = pubchem$AVERAGE_MASS, MONOISOTOPIC_MASS = pubchem$MONOISOTOPIC_MASS
           )


chemical_list_no_SMILES <- chemical_list_no_SMILES %>% filter(cas_number != chemicals_update$cas_number)


# recombine lists
chemical_list <- rbind(chemical_list_SMILES, chemicals_update)


# add comments and exclude those entries with remaining gaps
chemical_list <- chemical_list %>% mutate(EXCLUDE = ifelse(FOUND_BY == "No data retrieved from PubChem", 1, EXCLUDE),
                                          REMARKS = ifelse(FOUND_BY == "No data retrieved from PubChem", "no data", REMARKS))


write_csv(chemical_list, "C:/TEMP/EcoToxDB/test/pubchem_updated.csv")
