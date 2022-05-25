{rm(list = ls())
#require(data.table)
#    require(sqldf)
#require(tidyverse)

}

source("./R/Create_project.R")

# For documentation
# REcoTox 0.2.0

set.seed(4711456)

# Declare the database folder including the current version of EcoToX Knowledgebase
# A comprehensive name is for example "EcoTox_Fish_EC50"


database_path <- "c:/Data/UFZ_DATA/UFZ_Cloud/Databases/Ecotox/current"

# Declare the project folder to store the files of your query
project_path <- "c:/Data/UFZ_DATA/UFZ_Cloud/Projekte/EcoToxDB/Algae_EcoTox_220310_XX50"

# Declare the ecotox group ("Species Group")



# create the project
#
project <- create_project(database_path, project_path,
                          initalise_database_project = FALSE, # create the basic project from current ASCII files in DB folder
                          initalise_project = TRUE, # initializes the project folder
                          load_default = TRUE) # loads the default project in the project folder in the memory


# initialise_database_project (TRUE/FALSE) - create a new basic database out of the EcoTox ASCII files
# needs to be run each time the ASCII files are updated
#
# initialse_project (TRUE / FALSE): create the project folder
#
# load_default (TRUE / FALSE): copy the default database project from database folder to project folder

# Step 1: Run the first data preparation step to create the initial project

project <- prepare_data(project = project,
                        load_initial_project = TRUE,
                        new_project_path = NA,
                        save_project = FALSE
                        )

# Reload the results of the first step
# load(file.path(project_path, "initial_project.RData"))

# Step 2: Filter the data on the specified criteria
# Declare the settings for dataset filerting
dosing_group = "water_concentration" # i.e. mg/L (only available group in this version)
duration_d = c("d", "dph", "dpf")
duration_h = c("h", "ht", "hph", "hpf", "hbf", "hv")
duration_m = "mi"
ecotoxgroup = "Algae" # c("Algae", "Crustacean", "Fish")
effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV") # Algae/Fish
#effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV", "ITX") # Crustacean
habitat = "Water" #c("Non-Soil","Water","Soil")
kingdoms = NA # vector of specific algae kingdoms: c("Chromista","Plantae","Monera")
measurements = NA # vector of specific measurements
min_h = 0
min_d = 0
max_h = 120
max_d = 5
min_m = 0
max_m = 7200
species_selection = "all" # c("all", "manual", "standard_test_species")


# Run the workflow
project <- process_data(project,
                        dosing_group = dosing_group,
                        duration_d = duration_h,
                        duration_h = duration_h,
                        duration_m = duration_m,
                        ecotoxgroup = ecotoxgroup,
                        effects = effects,
                        habitat = habitat,
                        kingdoms = kingdoms,
                        measurements = measurements,
                        max_d = max_d,
                        min_d = min_d,
                        max_h = max_h,
                        min_h = min_h,
                        max_m = max_m,
                        min_m = min_m,
                        remove_formulation = FALSE,
                        save_project_steps = FALSE,
                        species_selection = species_selection
)
# A list of endpoints and species lists are stored in the project folder for selection

#load(file = file.path(project_path, paste0(ecotoxgroup,"_state1.RData")))

# Step 3: Read the modified lists in and process the data including unit conversion
# A list of chemicals is stored to update missing information on mol weights for data conversion

project <- process_data(project, save_project_steps = FALSE
)

#load(file = file.path(project_path, paste0(ecotoxgroup,"_state2.RData")))


# Step 4:
# Read in the the modified mol table and finally process the unit conversion
# This step creates a file named for example "fish_chemical_list.csv"
# Edit this list to include newly added compounds (imputation of phys.-
# chem. propertis and metadata)

project <- process_data(project, save_project_steps = FALSE
)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state3.RData")))

# Step 5:Process the final results and estimate the solubility domain
# Optional: Update the basic chemical list in the database folder

project <- process_data(project, save_project_steps = TRUE, update_chemicals = FALSE)


#load(file = file.path(project_path,paste0(ecotoxgroup,"_state4.RData")))
# save(project,file = file.path(project_path,paste0(prefix,"_pre_exclusion_project.RData")), compress = TRUE)


# Step 6: calculate and export the pivot table aggregating the results
project <- calculate_pivot_table(project = project, quantile = 0.05)



# do some final stuff
save(project, file = file.path(project_path,paste0(ecotoxgroup, "_processed_project.RData")), compress = TRUE)
rstudioapi::documentSave()
r_file <- rstudioapi::getSourceEditorContext()$path
file.copy(r_file,file.path(project_path), overwrite = TRUE)
