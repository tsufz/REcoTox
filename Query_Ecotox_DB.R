{rm(list = ls())
require(data.table)
#    require(sqldf)
require(tidyverse)

}

source("./R/Create_project.R")



set.seed(4711456)

# Declare the database folder including the current version of EcoToX Knowledgebase
# A comprehensive name is for example "EcoTox_Fish_EC50"


database_path <- "c:/Data/UFZ_DATA/UFZ_Cloud/Databases/Ecotox/current"

# Declare the project folder to store the files of your query
project_path <- "c:/Data/UFZ_DATA/UFZ_Cloud/Projekte/EcoToxDB/EcoToxDB_Fish_EC10_EC90"

# create the project
project <- create_project(database_path, project_path,
                          initalise_database_project = FALSE, # create the basic project from current ASCII files
                          initalise_project = TRUE, # initializes the project folder
                          load_default = TRUE) # loads the default project in the project folder

# initialise_database_project (TRUE/FALSE) - create a new basic database out of the EcoTox ASCII files
# needs to be run each time the ASCII files are updated
#
# initialse_project (TRUE / FALSE): create the project folder
#
# load_default (TRUE / FALSE): copy the default database project from database folder to project folder


# Declare effects of interest for the different species groups
# See EcoTox Knowledgebase documentation for the details
#

# Step 1: Run the first data preparation step to create the initial project

project <- prepare_data(project = project,
                        habitat = c("Water", "Non-Soil"),
                        effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV"), # Fish / Crustacean
                        # effects = c("MOR", "GRO", "POP", "REP"), # Algae
                        save_project = TRUE,
                        remove_formulation = FALSE,
                        new_project_path = NA,
                        load_initial_project = TRUE)

# Reload the results of the first step
load(file.path(project_path,"initial_project.RData"))



# Declare species group specific data

# Fish / Crustacean
#measurements = c("MORT", "SURV")

# Algae
#measurements = c("ABND","APCY","BMAS","CHLC","CHLO","CHLA","DBMS","DWGT","GPOP","GMOR","GGRO","INDX","MORT","PSYN","PSII","PPYT","PGRT","SPGR","SURV","VOLU","WGHT","WWGT")

# Declare the ecotox group ("Species Group")
#ecotoxgroup = "Algae"
#ecotoxgroup = "Fish"
ecotoxgroup = "Fish"
#ecotoxgroup = "Crustaceans"
#


# Step 2: Filter the data on the specified criteria
# A list of endpoints and species lists are stored in the project folder for selection

project <- process_data(project,
                        ecotoxgroup = ecotoxgroup,
                        max_h = 120,
                        max_d = 5,
                        kingdoms = c("Chromista","Plantae","Monera"), # valid only for algae
                        #species_selection = "selected")
                        #species_selection = "standard_species") # Standardized species only
                        all_species = TRUE, update_chemicals = TRUE
)


#load(file = file.path(project_path, paste0(ecotoxgroup,"_state1.RData")))

# Step 3: Read the modified lists in and process the data including unit conversion
# A list of chemicals is stored to update missing information on mol weights for data conversion

project <- process_data(project,
                        ecotoxgroup = ecotoxgroup,
                        max_h = 120,
                        max_d = 5,
                        kingdoms = c("Chromista","Plantae","Monera"), # valid only for algae
                        #species_selection = "selected")
                        #species_selection = "standard_species") # Standardized species only
                        all_species = TRUE, update_chemicals = TRUE
)


#load(file = file.path(project_path, paste0(ecotoxgroup,"_state2.RData")))


# Step 4:
# Read in the the modified mol table and finally process the unit conversion
# This step creates a file named for example "fish_chemical_list.csv"
# Edit this list to include newly added compounds (imputation of phys.-
# chem. propertis and metadata)

project <- process_data(project,
                        ecotoxgroup = ecotoxgroup,
                        max_h = 120,
                        max_d = 5,
                        kingdoms = c("Chromista","Plantae","Monera"), # valid only for algae
                        #species_selection = "selected")
                        #species_selection = "standard_species") # Standardized species only
                        all_species = TRUE, update_chemicals = TRUE
)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state3.RData")))

# Step 5:Process the final results and estimate the solubility domain
# Optional: Update the basic chemical list in the database folder

project <- process_data(project,
                        ecotoxgroup = ecotoxgroup,
                        max_h = 120,
                        max_d = 5,
                        kingdoms = c("Chromista","Plantae","Monera"), # valid only for algae
                        #species_selection = "selected")
                        #species_selection = "standard_species") # Standardized species only
                        all_species = TRUE,
                        update_chemicals = FALSE
)


#load(file = file.path(project_path,paste0(ecotoxgroup,"_state4.RData")))
# save(project,file = file.path(project_path,paste0(prefix,"_pre_exclusion_project.RData")), compress = TRUE)


# Step 6: calculate and export the pivot table aggregating the results
project <- calculate_pivot_table(project = project, quantile = 0.05)



# do some final stuff
save(project,file = file.path(project_path,paste0(ecotoxgroup,"_processed_project.RData")), compress = TRUE)
rstudioapi::documentSave()
r_file <- rstudioapi::getSourceEditorContext()$path
file.copy(r_file,file.path(project_path), overwrite = TRUE)
