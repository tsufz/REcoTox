library(tidyverse)
library(data.table)
project_path <- "c:/Data/UFZ_DATA/UFZ_Cloud/Projekte/EcoToxDB/EcoToxDB_Algae_NOEC_EC100"
files <- list.files(project_path, full.names = TRUE)
files
#fish <- read_csv(files[7]) # fish_chemical_list
chemical_list <- read_csv(files[1]) # fish_chemical_list
chemical_list_update <- read_csv(files[2])

#kek_JC <- read_csv(files[15])
#kek_OP <- read_csv(files[5])

kek_JC_base <- kek_JC %>% select(ID, Original_SMILES, Canonical_QSARr)
kek_JC_prop <- kek_JC %>% select(ID, LogP, LogD_74, logS, logS_74) %>%
    rename(JC_LOG_P = LogP, JC_LOG_D_74 = LogD_74, JC_LOG_S_74 = logS_74, JC_LOG_S = logS)

kek_OP_select <- kek_OP %>% select(ID, LogP_pred, AD_LogP, LogD74_pred, AD_LogD, LogWS_pred, AD_WS) %>%
    rename(OPERA_LOG_P = LogP_pred, OPERA_LOG_P_AD = AD_LogP, OPERA_LOG_D_74 = LogD74_pred,
           OPERA_LOG_D_AD = AD_LogD, OPERA_LOG_S_74 = LogWS_pred, OPERA_LOG_S_AD = AD_WS)

fish_merged <- fish %>% left_join(kek_JC_base, by = "ID")

fish_merged <- fish_merged %>% left_join(kek_OP_select, by = "ID")

fish_merged <- fish_merged %>% left_join(kek_JC_prop, by = "ID")

fish_ACD <- fish_merged %>% filter(!is.na(Canonical_QSARr)) %>% select(ID, Canonical_QSARr)

#write_csv(x = fish_ACD[1:1000, ], file.path(project_path, "ACD1_1000.csv"))
#write_csv(x = fish_ACD[1001:2000, ], file.path(project_path, "ACD1001_2000.csv"))
#write_csv(x = fish_ACD[2001:nrow(fish_ACD), ], file.path(project_path, "ACD2001_2829.csv"))


acd1 <- read_csv(file.path(project_path, "ACD1_1000.csv"))
acd2 <- read_csv(file.path(project_path, "ACD1001_2000.csv"))
acd3 <- read_csv(file.path(project_path, "ACD2001_2829.csv"))

acd <- rbind(acd1, acd2, acd3)

acd <- acd %>% select(ID, LogP, `LogD (pH = 7.40)`, `LogS (pH = 7.40)`) %>% rename(ID = ID, ACD_LOG_P = LogP,
                                                                                   ACD_LOG_D_74 = 'LogD (pH = 7.40)',
                                                                                   ACD_LOG_S_74 = 'LogS (pH = 7.40)')

fish_merged <- fish_merged %>% left_join(acd, by = "ID")
write_csv(fish_merged, file = file.path(project_path, "fish_chemicals_merged.csv"))

fish_chemicals <- read_csv(file.path(project_path, "fish_chemical_list.csv"))

chemprop <- data.table(fish)
chemical_list <- data.table(fish_merged)
#chemical_list <- format_chemical_properties(chemical_list)
chemprop <- data.table(merge(chemprop, chemical_list, by = 'cas_number', all = TRUE, suffixes = c("", ".update")))
col_names1 <- grep(".update", colnames(chemprop))
col_names2 <- colnames(chemprop)[grep(".update", colnames(chemprop))]
col_names3 <- gsub(".update", replacement = "", col_names2)



length_progressbar <- length(col_names1)
pb <- progress::progress_bar$new(
    format = "[EcoToxR]:  Updating the chemical properties [:bar] :percent ETA: :eta",
    total = length_progressbar, clear = FALSE, width = 80)


for(i in 1:length(col_names1)){
    c1 <- which(colnames(chemprop) == col_names3[i])
    c2 <- which(colnames(chemprop) == col_names2[i])
    for(j in 1:nrow(chemprop)){
        if(!is.na(chemprop[j,..c2])){
            chemprop[j, c1] <- chemprop[j, ..c2]
        } else next()
    }
    pb$tick()
}

chemprop <- data.table(chemprop)
chemprop[,(col_names2) := NULL]

write_csv(x = chemprop, file = file.path(project_path, "check_new_chemicals.csv"))



files <- list.files(project_path, full.names = TRUE)
files
#fish <- read_csv(files[7]) # fish_chemical_list
chemprop <- read_csv(files[1]) # fish_chemical_list
chemical_list <- read_csv(files[2])

chemical_list <- chemical_list %>% left_join(chemprop %>% select(MoleculeID, cas_number))


chemprop <- data.table(chemprop)
chemical_list <- data.table(chemical_list)
#chemical_list <- format_chemical_properties(chemical_list)
chemprop <- data.table(merge(chemprop, chemical_list, by = 'cas_number', all = TRUE, suffixes = c("", ".update")))
col_names1 <- grep(".update", colnames(chemprop))
col_names2 <- colnames(chemprop)[grep(".update", colnames(chemprop))]
col_names3 <- gsub(".update", replacement = "", col_names2)



length_progressbar <- length(col_names1)
pb <- progress::progress_bar$new(
    format = "[EcoToxR]:  Updating the chemical properties [:bar] :percent ETA: :eta",
    total = length_progressbar, clear = FALSE, width = 80)


for(i in 1:length(col_names1)){
    c1 <- which(colnames(chemprop) == col_names3[i])
    c2 <- which(colnames(chemprop) == col_names2[i])
    for(j in 1:nrow(chemprop)){
        if(!is.na(chemprop[j,..c2])){
            chemprop[j, c1] <- chemprop[j, ..c2]
        } else next()
    }
    pb$tick()
}

chemprop <- data.table(chemprop)
chemprop[,(col_names2) := NULL]

write_csv(x = chemprop, file = file.path(project_path, "check_new_chemicals.csv"))

