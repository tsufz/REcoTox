#' @importFrom dplyr mutate left_join filter %>% select group_by ungroup
#' @importFrom dplyr pull arrange bind_rows rowwise if_else
#' @importFrom dplyr case_when slice_min mutate_all na_if
#' @importFrom tibble tibble add_column add_row
#' @importFrom readr read_delim write_csv read_csv
#' @importFrom data.table %like%
#' @importFrom purrr is_empty
#' @importFrom progress progress_bar
#' @importFrom webchem pc_prop
#'
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\nThis is REcoTox version", packageVersion("REcoTox"), "\n"))
}
