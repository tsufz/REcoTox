{rm(list = ls())
#require(data.table)
#    require(sqldf)
#require(tidyverse)

}
# For documentation
# REcoTox 0.3.5

# Use only for development purposes
# source("./R/Create_project.R")

library(REcoTox)

# Declare the database folder including the current version of EcoToX Knowledgebase
# A comprehensive name is for example "EcoTox_Fish_EC50"


database_path <- "path_to_ecotox_unzipped_ascii_files"

# Declare the project folder to store the files of your query
project_path <- "path_to_project_folder"

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
# The filtering needs some knowledge on the internal structures of the database
# For reference see https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix
dosing_group = "water_concentration" # i.e. mg/L (only available group in this version)

# time based settings
duration_d = c("d", "dph", "dpf") # day based units
duration_h = c("h", "ht", "hph", "hpf", "hbf", "hv") # hour based units
duration_m = "mi" # minute based units
min_h = 0 # minimum hours
min_d = 0 # minimum days
min_m = 0 # minimum minutes

max_h = 120 # maximum hours
max_d = 5 # maximum days
max_m = 7200 # maximum minutes

# species base settings
ecotoxgroup = "Algae" # c("Algae", "Crustacean", "Fish")
species_selection = "all" # c("all", "manual", "standard_test_species")
habitat = "Water" #c("Non-Soil","Water","Soil")
kingdoms = NA # vector of specific algae kingdoms: c("Chromista","Plantae","Monera")

# effects and measurements
effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV") # Algae/Fish
#effects = c("MOR", "GRO", "POP", "REP", "MPH", "DEV", "ITX") # Crustacean

measurements = NA # vector of specific measurements for refined selection

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
# Optional: Save the project file to the project folder.

project <- process_data(project, save_project_steps = FALSE
)

#load(file = file.path(project_path,paste0(ecotoxgroup,"_state3.RData")))

# Step 5:Process the final results and estimate the solubility domain
# Optional: Update the basic chemical list in the database folder
# The update is recommended, if the chemical list was edited / updated
# Optional: Save the project to the project_folder

project <- process_data(project, save_project_steps = FALSE, update_chemicals = FALSE)


#load(file = file.path(project_path,paste0(ecotoxgroup,"_state4.RData")))
# save(project,file = file.path(project_path,paste0(prefix,"_pre_exclusion_project.RData")), compress = TRUE)


# Step 6: calculate and export the pivot table aggregating the results
# Select the value for the percentile cut-off of the ECx values
project <- calculate_pivot_table(project = project, quantile = 0.05)



# do some final stuff
save(project, file = file.path(project_path,paste0(ecotoxgroup, "_processed_project.RData")), compress = TRUE)
rstudioapi::documentSave() # save current R script
r_file <- rstudioapi::getSourceEditorContext()$path
file.copy(r_file,file.path(project_path), overwrite = FALSE) # copy the script to the project_folder
