geo_acquisition_method_for_role <- function(role) {
  if (identical(role, "official_quick_metadata")) return("GEOquery::getGEOfile")
  if (identical(role, "series_matrix")) return("GEOquery::getGEO")
  if (identical(role, "official_supplement_filelist")) return("curl_auxiliary_filelist")
  if (identical(role, "raw_archive")) return("retain_existing_source_no_auto_download")
  "GEOquery::getGEOSuppFiles"
}

is_analysis_matrix_role <- function(role) {
  role %in% c("series_matrix", "non_normalized_expression", "raw_counts")
}

lightweight_availability <- function(cohort, cohort_dir) {
  items <- cohort$downloads
  paths <- vapply(items, function(x) file.path(cohort_dir, x$filename), character(1))
  roles <- vapply(items, `[[`, character(1), "role")
  methods <- vapply(roles, geo_acquisition_method_for_role, character(1))
  analysis_indices <- vapply(roles, is_analysis_matrix_role, logical(1))
  required_indices <- methods != "retain_existing_source_no_auto_download"
  list(
    analysis_matrix_present = any(file.exists(paths[analysis_indices])),
    missing_required_non_source = paths[required_indices & !file.exists(paths)],
    missing_retained_source = paths[!required_indices & !file.exists(paths)]
  )
}
