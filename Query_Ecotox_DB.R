{rm(list = ls())
require(data.table)
    require(sqldf)
    require(tidyverse)

}

source("./R/EcoToxDB/Create_project.R")



set.seed(4711456)
#source('C:/git/Environmental_Priority_Mixtures/R/helperfunctions.R')
#database_path <- "Y:/Home/schulzet/UFZ/Databases/Ecotox/current" #MSG
database_path <- "d:/UFZ_DATA/UFZ/Databases/Ecotox/Current" #WANA61
#database_path <- "d:/Daten/UFZ/UFZ_Data/Databases/Ecotox/current" #home

#project_path <- "Y:/Home/schulzet/UFZ/Projekte/EcoToxDB/EcoTox_Algae_selected_species_XX50"
project_path <- "d:/UFZ_DATA/UFZ/Projekte/EcoToxDB/Chironimus_EC50"
#project_path <- "d:/daten/UFZ/UFZ_Data/Projekte/EcoToxDB/EcoTox_Algae_selected_species_XX50"


# create the project



# create the project
project <- create_project(database_path, project_path,
                          initalise_database_project = FALSE,
                          initalise_project = TRUE,
                          load_default = TRUE)


#effects <- c("GRO","DEV","MPH","MOR","POP") # Fish
#effects <- c("GRO","MPH","MOR","POP") # Algae, Crustacean
#effect = c("MOR","GRO","DEV")
#habitat = c("Non-Soil","Water","Soil")

project <- prepare_data(project, habitat = c("Water","Non-Soil"),
                        #effects = c("MOR","GRO","POP"),
                        effects = c("GRO","DEV","MPH","MOR","POP"),
                        remove_formulation = TRUE)

load(file.path(project_path,"initial_project.RData"))
# Process group specific data

#measurement = c("MORT","GGRO","SURV","GMOR")
# Fish
measurements = c("MORT","SURV")
# Crustacean / Insects

# Algae
#measurements = c("ABND","APCY","BMAS","CHLC","CHLO","CHLA","DBMS","DWGT","GPOP","GMOR","GGRO","INDX","MORT","PSYN","PSII","PPYT","PGRT","SPGR","SURV","VOLU","WGHT","WWGT")

#ecotoxgroup = "Algae"
#ecotoxgroup = "Fish"
#ecotoxgroup = "Crustacean"
ecotoxgroup = "Insects" # Chironimus included here
#ecotoxgroup = "Worms" # Nematods
#ecotoxgroup = "Mammals"
#ecotoxgroup = "Reptiles"
#ecotoxgroup = "Amphibians"
#ecotoxgroup = "Ivertebrates" # Other invertebrates



# Prepare the endpoint and species lists for edits

project <- process_data(project, ecotoxgroup = ecotoxgroup,
                        max_h = 120, max_d = 5,
                        measurements = measurements)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state1.RData")))

# Read the modified lists in and process the data including unit conversion

project <- process_data(project, ecotoxgroup = ecotoxgroup,
                        max_h = 120, max_d = 5,
                        measurements = measurements)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state2.RData")))

# Read in the the modified mol table and finally process the unit conversion

project <- process_data(project, ecotoxgroup = ecotoxgroup,
                        max_h = 120, max_d = 5,
                        measurements = measurements)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state3.RData")))

# Process the final results and estimate the solubility domain
project <- process_data(project, ecotoxgroup = ecotoxgroup,
                        max_h = 120, max_d = 5,
                        measurements = measurements)

#project2 -> project

#project$object$state
#object <- project$object

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state3.RData")))
#save(project,file = file.path(project_path,paste0(prefix,"_pre_exclusion_project.RData")), compress = TRUE)


# calculate and export the pivot table
project <- calculate_pivot_table(project = project, quantile = 0.05)



# do some final stuff
save(project,file = file.path(project_path,paste0(ecotoxgroup,"_processed_project.RData")), compress = TRUE)
#load(file.path(project_path,"Fish_selected_pre_exclusion_project.RData"))


rstudioapi::documentSave()
r_file <- rstudioapi::getSourceEditorContext()$path
file.copy(r_file,file.path(project_path),overwrite = TRUE)
