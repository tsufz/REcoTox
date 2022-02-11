source("./R/helperfunctions.R")
source("./R/Prepare_data.R")
source("./R/Calculate_pivot.R")
#source("./R/EcoToxDB/Create_project.R")
source("./R/Process_data.R")
require(progress)

create_project <- function(database_path, project_path, initalise_database_project = FALSE,
                           initalise_project = FALSE, load_default = FALSE){

  # initialise a new basic project, otherwise just load and update the project file
  if (initalise_database_project == TRUE) {
    project <- list()
    project$object <- list()
    project$database_path <- suppressWarnings(normalizePath(database_path))
    object <- project$object
    project$files <- normalizePath(list.files(database_path, pattern = "*.txt", full.names = T, recursive = T))
    message("[EcoToxR]:  Initializing a new project in the database folder. The old will be overwritten.")

    # read the file list and harmonise use of NA
    message("[EcoToxR]:  Read the Ecotox Knowledgebase ASCII files.")
    suppressWarnings({
      message("[EcoToxR]:  Read tests")
      object$tests <- read_delim(project$files[grep("tests.txt", project$files)], delim = "|", na = c("NA","","NR","--","NC","/"), quote = "\"")

      message("[EcoToxR]:  Read chemicals.")
      object$chemicals <- read_delim(project$files[grep("chemicals.txt", project$files)], delim = "|", na = c("NA","","NR","--","NC","/"), quote = "\"")

      object$chemprop <- dplyr::tibble(create_chemical_properties(project$database_path))

      message("[EcoToxR]:  Read species.")
      object$species <- read_delim(project$files[grep("species.txt", project$files)], delim = "|", na = c("NA","","NR","--","NC","/"), quote = "\"")

      message("[EcoToxR]:  Read results.")
      #object$results <- fread(project$files[grep("results.txt", project$files)], sep = "|", header = T, na.strings = c("NA","","NR","--","NC","/"), dec = ".", stringsAsFactors = FALSE)
      object$results <- read_delim(project$files[grep("results.txt", project$files)], delim = "|", na = c("NA","","NR","--","NC","/"), quote = "\"")

      message("[EcoToxR]:  Read references.")
      object$references <- read_delim(project$files[grep("references.txt", project$files)], delim = "|", na = c("NA","","NR","--","NC","/"), quote = "\"")
    })

    #  Trim lists
    names(object$chemicals)[which(names(object$chemicals) %like% "ecotox_group")] <- "compound_class"
    names(object$chemicals)[which(names(object$chemicals) %like% "dtxsid")] <- "DTXSID"

    # harmonize name for query with chemicals
    names(object$tests)[which(names(object$tests) %like% "test_cas")] <- "cas_number"

    # Add a column with a CAS in ususal format
    object$chemicals <- object$chemicals %>% mutate(cas = paste0(substr(as.character(object$chemicals$cas_number), 1, nchar(as.character(object$chemicals$cas_number)) - 3), "-",
                                                                 substr(as.character(object$chemicals$cas_number), nchar(as.character(object$chemicals$cas_number)) - 2, nchar(as.character(object$chemicals$cas_number)) - 1), "-",
                                                                 substr(as.character(object$chemicals$cas_number), nchar(as.character(object$chemicals$cas_number)), nchar(as.character(object$chemicals$cas_number)))))

    message("[EcoToxR]:  Saving the basic project in the database folder.")
    project$object <- object
    object <- NULL
    save(project, file = file.path(database_path,"project.RData"), compress = TRUE)
  }


  if (initalise_project == TRUE) {
    # Check the directory
    if (!dir.exists(project_path)) {
      message("[EcoToxR]:  Creating the project directory.")
      dir.create(project_path)
    } else if (!is_empty(project_path)) {
      message("[EcoToxR]:  The project directory is not empty. Please checkout content")
      message("[EcoToxR]:  and delete it or use another path before you continue.")
    } else if (is_empty(project_path)) {
      message("[EcoToxR]:  The project directory is empty.")
    }
    if (!exists("project") & load_default == FALSE) {
      load_default = TRUE
    }
    if (exists("project") & load_default == FALSE) {
      project$object$chemprop <- dplyr::tibble(create_chemical_properties(project$database_path))
      project$project_path <- suppressWarnings(normalizePath(project_path))
    }
  }

  if (load_default == TRUE) {
    message("[EcoToxR]:  Loading the default project file.")
    load(file.path(database_path,"project.RData"))
    project$database_path <- suppressWarnings(normalizePath(database_path))
    project$project_path <- suppressWarnings(normalizePath(project_path))
    project$object$chemprop <- create_chemical_properties(database_path)
    project$files <- suppressWarnings(normalizePath(list.files(database_path, pattern = "*.txt", full.names = T, recursive = T)))
  }

  return(project)
}
