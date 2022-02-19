library(webchem)
library(progress)
library(tidyverse)
library(data.table)

chemical_list <- read_csv("c:/TEMP/cheos/opera/180222_opera_input_file.csv")
opera_predicted <- read_csv("c:/temp/cheos/opera/ECOTOXDB_151222_opera_1000_OPERA_2.7_prd.csv")

opera_predicted_select <- opera_predicted %>% select(MoleculeID, LogP_pred, AD_LogP, LogD74_pred, AD_LogD, LogWS_pred, AD_WS) %>%
    rename(c("Molecule_ID_OPERA" = "MoleculeID", "OPERA_LOG_P" = "LogP_pred", "AD_OPERA_LOG_P" = "AD_LogP",
             "OPERA_LOG_D_74" = "LogD74_pred", "AD_OPERA_LOG_D_74" = "AD_LogD", "OPERA_LOG_S" = "LogWS_pred",
             "AD_OPERA_LOG_S" = "AD_WS"))


chemical_list_merge <- chemical_list %>% left_join(opera_predicted_select, by = "Molecule_ID_OPERA")
chemical_list_acd <- chemical_list_merge %>% filter(is.na(EXCLUDE)) %>% select(MoleculeID, cas_number, QSAR_READY_SMILES) %>% arrange(cas_number)

acd_1000 <- chemical_list_acd %>% slice(n = 1:1000)
acd_2000 <- chemical_list_acd %>% slice(n = 1001:2000)
acd_3000 <- chemical_list_acd %>% slice(n = 2001:3000)
acd_4000 <- chemical_list_acd %>% slice(n = 3001:4000)


write_csv(acd_1000, "c:/TEMP/cheos/opera/180222_acd_1000.csv")
write_csv(acd_2000, "c:/TEMP/cheos/opera/180222_acd_2000.csv")
write_csv(acd_3000, "c:/TEMP/cheos/opera/180222_acd_3000.csv")
write_csv(acd_4000, "c:/TEMP/cheos/opera/180222_acd_4000.csv")

acd_1000 <- read_csv("c:/TEMP/cheos/opera/180222_acd_1000.csv")
acd_2000 <- read_csv("c:/TEMP/cheos/opera/180222_acd_2000.csv")
acd_3000 <- read_csv("c:/TEMP/cheos/opera/180222_acd_3000.csv")
acd_4000 <- read_csv("c:/TEMP/cheos/opera/180222_acd_4000.csv")

acd <- rbind(acd_1000, acd_2000, acd_3000, acd_4000)

acd <- acd %>% select(MoleculeID, LogP, `LogD (pH = 7.40)`, `LogS (pH = 7.40)`) %>% rename("ACD_LOG_P" = "LogP",
                                                                                                       "ACD_LOG_D_74" = "LogD (pH = 7.40)",
                                                                                                       "ACD_LOG_S_74" = "LogS (pH = 7.40)")
chemical_list_final <- chemical_list_merge %>% left_join(acd, by ="MoleculeID")


chemicals_final <- chemical_list_final %>% select(cas_number, QSAR_READY_SMILES, OPERA_LOG_P, AD_OPERA_LOG_P, OPERA_LOG_D_74, AD_OPERA_LOG_D_74,
                                                  OPERA_LOG_S, AD_OPERA_LOG_S, ACD_LOG_P, ACD_LOG_D_74, ACD_LOG_S_74, JC_LOG_P, JC_LOG_D_74,
                                                  JC_LOG_S_74, EXCLUDE, REMARKS) %>% arrange(across(cas_number))

write_csv(chemicals_final, "c:/TEMP/cheos/opera/180222_chemicals_properties.csv")
