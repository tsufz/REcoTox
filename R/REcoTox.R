#' REcoTox - a workflow to process US EPA ECOTOX Knowledgebase ASCII files
#'
#' The search and extraction of experimental ecotoxicological information is
#' often a tedious work. A good and comprehensive data source is the US EPA
#' ECOTOX Knowledgebase \insertCite{Olker_2022}{REcoTox}.
#' It contains about 1.1 million data points for more than 12,000 chemicals
#' and 13,000 single species. However, for a high-throughput hazard assessment,
#' it is not possible to extract all relevant data of the online database.
#'
#' The purpose of REcoTox is to extract the relevant information and to
#' aggregate the data based on the user criteria out of the entire
#' database ASCII files.
#'
#' @section REcoTox functions:
#' The REcoTox functions:
#' \itemize{
#' \item \code{\link{create_project}}
#' Initialise the default project and read in EcoTox Knowlegdebase \code{ascii}
#' files from the \code{database_folder}.
#' \item \code{\link{prepare_data}}
#' Prepares the initial project and merges the tables for downstream queries.
#' \item \code{\link{process_data}}
#' Processing and filtering workflow.
#' \item \code{\link{aggregate_results}}
#' Aggregates data and creates the final result table.
#' }
#'
#' @docType package
#' @name REcoTox
#'
#' @references
#' \insertRef{Olker_2022}{REcoTox}
#'
NULL
