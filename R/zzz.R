#' @importFrom dplyr mutate left_join filter %>% select group_by ungroup
#' @importFrom dplyr pull add_column arrange row_bind rowwise if_else
#' @importFrom dplyr ungroup case_when add_row slice_min mutate_all na_if
#' @importFrom tibble as_tibble
#' @importFrom readr read_delim write_csv read_csv
#' @importFrom data.table %like%
#' @importFrom purrr is_empty
#' @importFrom progress progress_bar
#' @importFrom webchem pc_prop


# globalVariables(c(":=", "StdName", "Class", "RT", "MZmineID", "Name",
#                 "Concentration", "project_folder", "MZ", "FREQ",
#                  "Filename", "IS_gap_filling", "MZ_ppm_error",
#                  "Row_MZ", "Row_RT", "Substance_Name", "blank_peaks_detected",
#                  "blank_peaks_max", "peak_method", "blank_peaks_mean",
#                  "blank_peaks_min", "blank_peaks_sd", "blank_peaks_threshold",
#                  "calibration_aligned_data", "calibration_data_by_concentration",
#                  "calibration_peaks_detected", "compound_IS_assignment",
#                  "deltaMZ", "deltaRT", "in_cal_num", "localmax_pos",
#                  "localmin_pos", "missing_targets", "new_Class",
#                  "new_StdName", "peaklist_calibration_clean",
#                  "peaklist_samples_clean", "project_folder",
#                  "remarks", "remove_annotation", "sample_peaks_detected",
#                  "total_peaks_detected", "calibration_monocity",
#                  "calibration_column_names", "qm_prediction_func_all", "x"))
