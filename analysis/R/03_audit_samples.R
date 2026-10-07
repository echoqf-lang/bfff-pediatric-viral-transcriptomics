options(stringsAsFactors = FALSE)

required_root_files <- c(
  file.path("analysis", "config", "cohorts.yml"),
  file.path("data", "metadata", "all_samples_raw.tsv")
)
if (!all(file.exists(required_root_files))) {
  stop("Run Task 4 from the project root; required frozen metadata/configuration is missing.")
}

trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

after_colon <- function(x) {
  x <- trim_na(x)
  out <- sub("^[^:]+:[[:space:]]*", "", x)
  out[is.na(x)] <- NA_character_
  out
}

numeric_after_colon <- function(x) suppressWarnings(as.numeric(after_colon(x)))

normalize_sex <- function(x) {
  value <- tolower(trim_na(after_colon(x)))
  out <- ifelse(value %in% c("m", "male"), "male",
                ifelse(value %in% c("f", "female"), "female", NA_character_))
  out
}

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "",
    fileEncoding = "UTF-8"
  )
}

read_series_matrix_header <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  accession <- platform <- character()
  sample_ids <- character()
  stopped <- FALSE
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (identical(line, "!series_matrix_table_begin")) {
      stopped <- TRUE
      break
    }
    if (startsWith(line, "!Series_geo_accession")) {
      accession <- gsub('"', "", trimws(sub("^[^\t]+\t", "", line)))
    } else if (startsWith(line, "!Series_platform_id")) {
      platform <- c(platform, gsub('"', "", trimws(sub("^[^\t]+\t", "", line))))
    } else if (startsWith(line, "!Series_sample_id")) {
      value <- gsub('"', "", trimws(sub("^[^\t]+\t", "", line)))
      sample_ids <- strsplit(value, "[[:space:]]+")[[1L]]
      sample_ids <- sample_ids[nzchar(sample_ids)]
    }
  }
  list(
    accession = accession,
    platforms = unique(platform),
    sample_ids = sample_ids,
    stopped_before_table = stopped
  )
}

read_quick_soft_accession <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hits <- lines[startsWith(lines, "!Series_geo_accession =")]
  if (length(hits) != 1L) stop("Quick SOFT must contain exactly one Series accession: ", path)
  trimws(sub("^!Series_geo_accession[[:space:]]*=[[:space:]]*", "", hits))
}

blank_manifest <- function(raw) {
  data.frame(
    series_accession = raw$series_accession,
    sample_id = raw$geo_accession,
    title = raw$title,
    analysis_role = raw$analysis_role,
    outcome_previously_observed = as.logical(raw$outcome_previously_observed),
    subject_id = paste(raw$series_accession, raw$geo_accession, sep = "_SUBJ_"),
    center = NA_character_, platform = raw$platform_id,
    analysis_unit = NA_character_, timepoint = NA_character_, status = NA_character_,
    pathogen = NA_character_, age_months = NA_real_, sex = NA_character_,
    severity = NA_character_, coinfection = NA_character_, treatment = NA_character_,
    batch = NA_character_, tissue = NA_character_,
    replicate_type = "none", replicate_group = NA_character_, include = FALSE,
    inclusion_scope = "excluded", exclusion_reason = "unclassified_conservative_exclusion",
    needs_manual_review = TRUE, metadata_evidence = NA_character_,
    metadata_discrepancy = NA_character_,
    meta_eligible_confirmatory = FALSE,
    severity_analysis_eligible = FALSE,
    processed_only_qc_boundary = "metadata_frozen; CEL/FASTQ-level QC not claimed",
    stringsAsFactors = FALSE
  )
}

audit_samples <- function(raw) {
  stopifnot(!anyDuplicated(raw$geo_accession))
  m <- blank_manifest(raw)

  set_rows <- function(idx, ...) {
    values <- list(...)
    for (name in names(values)) m[idx, name] <<- values[[name]]
  }

  # GSE105450: the public series design declares five healthy technical repeats.
  # They are resolved only by identical sex/race/date/age/condition metadata across
  # Discovery and Validation and never by expression similarity.
  idx <- raw$series_accession == "GSE105450"
  condition <- after_colon(raw$characteristics_ch1__repeat_4[idx])
  severity <- ifelse(condition == "RSV Inpatient", "hospitalized",
                     ifelse(condition == "RSV Outpatient", "outpatient", "healthy"))
  status <- ifelse(condition == "Healthy", "healthy_control", "RSV_case")
  local_rows <- which(idx)
  set_rows(
    idx, center = "Columbus_NCH", timepoint = "acute_primary", status = status,
    pathogen = ifelse(status == "RSV_case", "RSV", "none"),
    age_months = numeric_after_colon(raw$characteristics_ch1__repeat_3[idx]),
    sex = normalize_sex(raw$characteristics_ch1[idx]), severity = severity,
    coinfection = "none", batch = paste(
      after_colon(raw$characteristics_ch1__repeat_5[idx]),
      after_colon(raw$characteristics_ch1__repeat_6[idx]), sep = ";"
    ), tissue = "whole_blood", include = TRUE,
    inclusion_scope = "confirmatory_primary", exclusion_reason = NA_character_,
    needs_manual_review = FALSE,
    metadata_evidence = "condition+tissue+age+sex+race+collection_date; series design confirms acute RSV/healthy and five technical repeats"
  )
  healthy_local <- status == "healthy_control"
  duplicate_key <- paste(
    raw$characteristics_ch1[idx], raw$characteristics_ch1__repeat_1[idx],
    raw$characteristics_ch1__repeat_2[idx], raw$characteristics_ch1__repeat_3[idx],
    raw$characteristics_ch1__repeat_4[idx], sep = "||"
  )
  pair_sizes <- table(duplicate_key[healthy_local])
  technical_keys <- names(pair_sizes[pair_sizes == 2L])
  technical_local <- healthy_local & duplicate_key %in% technical_keys
  if (length(technical_keys) != 5L || sum(technical_local) != 10L) {
    stop("GSE105450 public metadata no longer resolves exactly five healthy technical-replicate pairs.")
  }
  for (key in technical_keys) {
    pair_local <- which(idx)[duplicate_key == key & healthy_local]
    pair_local <- pair_local[order(m$sample_id[pair_local], method = "radix")]
    group <- paste0("GSE105450_TECH_", m$sample_id[pair_local[[1L]]])
    m$subject_id[pair_local] <- paste0("GSE105450_SUBJ_", m$sample_id[pair_local[[1L]]])
    m$replicate_type[pair_local] <- "technical"
    m$replicate_group[pair_local] <- group
    m$include[pair_local[[2L]]] <- FALSE
    m$inclusion_scope[pair_local[[2L]]] <- "excluded"
    m$exclusion_reason[pair_local[[2L]]] <- "technical_replicate"
  }

  # GSE103842: all 62 acute hospitalized RSV and 12 asymptomatic controls.
  idx <- raw$series_accession == "GSE103842"
  condition <- after_colon(raw$characteristics_ch1__repeat_2[idx])
  status <- ifelse(condition == "Control", "healthy_control", "RSV_case")
  set_rows(
    idx, center = NA_character_, timepoint = "acute_primary", status = status,
    pathogen = ifelse(status == "RSV_case", "RSV", "none"),
    age_months = numeric_after_colon(raw$characteristics_ch1[idx]),
    sex = normalize_sex(raw$characteristics_ch1__repeat_1[idx]),
    severity = ifelse(status == "RSV_case", "hospitalized", "healthy"),
    coinfection = "none", batch = sub("_[A-L]$", "", raw$description[idx]),
    tissue = "whole_blood", include = TRUE,
    inclusion_scope = "confirmatory_primary", exclusion_reason = NA_character_,
    needs_manual_review = FALSE,
    metadata_evidence = "condition1+condition2+RSV genotype+tissue; series design states blood within 24h and asymptomatic controls"
  )

  # GSE188427: subject prefix links Day 1/30/180; only Day 1 can substitute.
  idx <- raw$series_accession == "GSE188427"
  title_prefix <- sub("_.*$", "", raw$title[idx])
  subject_prefix <- sub("-D(1|30|180)$", "", title_prefix)
  day <- after_colon(raw$characteristics_ch1__repeat_1[idx])
  arm <- after_colon(raw$characteristics_ch1__repeat_2[idx])
  source_arm <- sub(".*_", "", raw$source_name_ch1[idx])
  title_arm <- sub(".*_", "", raw$title[idx])
  arm_conflict <- source_arm != arm | title_arm != arm
  disease <- after_colon(raw$characteristics_ch1[idx])
  status <- ifelse(disease == "Healthy person", "healthy_control", "RSV_case")
  timepoint <- c("Day 1" = "day_1", "Day 30" = "day_30", "Day 180" = "day_180")[day]
  include <- timepoint == "day_1"
  set_rows(
    idx, subject_id = paste0("GSE188427_SUBJ_", subject_prefix), center = NA_character_,
    timepoint = unname(timepoint), status = status,
    pathogen = ifelse(status == "RSV_case", "RSV", "none"),
    severity = ifelse(
      arm_conflict, NA_character_,
      ifelse(arm == "IN", "hospitalized", ifelse(arm == "OUT", "outpatient", "healthy"))
    ),
    coinfection = "none", batch = raw$description[idx], tissue = "whole_blood",
    replicate_type = "none", replicate_group = NA_character_,
    include = include,
    inclusion_scope = ifelse(include, "substitution_sensitivity_only", "excluded"),
    exclusion_reason = ifelse(include, NA_character_, "non_primary_longitudinal_timepoint"),
    needs_manual_review = arm_conflict,
    metadata_evidence = "title subject prefix+day+disease; arm used for severity only when title, source_name, and characteristics agree",
    metadata_discrepancy = ifelse(
      source_arm != arm,
      "source_name arm IN conflicts with title and characteristics arm OUT; severity unresolved",
      NA_character_
    )
  )
  g188_rows <- which(idx)
  g188_counts <- table(m$subject_id[g188_rows])
  repeated <- unname(g188_counts[m$subject_id[g188_rows]]) > 1L
  m$replicate_type[g188_rows[repeated]] <- "longitudinal"
  m$replicate_group[g188_rows[repeated]] <- paste0("GSE188427_LONG_", subject_prefix[repeated])

  # GSE38900: retain only explicitly non-NCH Dallas/Turku acute RSV and healthy.
  idx <- raw$series_accession == "GSE38900"
  source <- raw$source_name_ch1[idx]
  center <- ifelse(grepl("Columbus", source), "Columbus_NCH",
                   ifelse(grepl("Finnish", source), "Turku",
                          ifelse(grepl("Dallas|Dalls", source), "Dallas", NA_character_)))
  diagnosis <- raw$characteristics_ch1__repeat_2[idx]
  status <- ifelse(grepl("healthy", diagnosis, ignore.case = TRUE), "healthy_control",
                   ifelse(grepl("follow|1-2 month", paste(source, diagnosis), ignore.case = TRUE), "recovery",
                          ifelse(grepl("RSV", diagnosis), "RSV_case", "other_virus_case")))
  timepoint <- ifelse(status == "recovery", "recovery_1_2_month", "acute_primary")
  normalized_subject <- gsub("[^0-9]", "", sub("RSV.*$", "", raw$description[idx]))
  normalized_subject[!nzchar(normalized_subject)] <- raw$geo_accession[idx][!nzchar(normalized_subject)]
  subject_key <- ifelse(
    !is.na(center), paste("GSE38900", center, normalized_subject, sep = ":"),
    paste("GSE38900", "unknown_center", raw$geo_accession[idx], sep = ":")
  )
  include <- center %in% c("Dallas", "Turku") & status %in% c("RSV_case", "healthy_control") & timepoint == "acute_primary"
  reason <- ifelse(timepoint != "acute_primary", "recovery_timepoint",
                   ifelse(status == "other_virus_case", "non_RSV_virus",
                          ifelse(center == "Columbus_NCH", "NCH_overlap_risk_excluded", "ambiguous_center_or_group")))
  reason[include] <- NA_character_
  set_rows(
    idx, subject_id = subject_key, center = center,
    analysis_unit = ifelse(!is.na(center), paste(raw$platform_id[idx], center, sep = "__"), NA_character_),
    timepoint = timepoint, status = status,
    pathogen = ifelse(status == "RSV_case" | status == "recovery", "RSV",
                      ifelse(status == "healthy_control", "none", after_colon(diagnosis))),
    age_months = numeric_after_colon(raw$characteristics_ch1[idx]),
    sex = normalize_sex(raw$characteristics_ch1__repeat_1[idx]),
    severity = ifelse(status %in% c("RSV_case", "other_virus_case", "recovery"), "hospitalized", "healthy"),
    coinfection = "none", tissue = "whole_blood",
    replicate_type = "none", replicate_group = NA_character_,
    include = include, inclusion_scope = ifelse(include, "exploratory_previously_observed", "excluded"),
    exclusion_reason = reason, needs_manual_review = is.na(center),
    metadata_evidence = "source_name explicitly identifies Columbus/Dallas/Finnish center and acute/follow-up/healthy group"
  )
  g389_rows <- which(idx)
  g389_counts <- table(m$subject_id[g389_rows])
  repeated <- unname(g389_counts[m$subject_id[g389_rows]]) > 1L
  m$replicate_type[g389_rows[repeated]] <- "longitudinal"
  m$replicate_group[g389_rows[repeated]] <- paste0("GSE38900_LONG_", m$subject_id[g389_rows[repeated]])

  # GSE155925: hospitalized cross-virus/severity exploration; negative and mixed-virus samples excluded.
  idx <- raw$series_accession == "GSE155925"
  pathogen <- after_colon(raw$characteristics_ch1__repeat_3[idx])
  mixed <- grepl(",", pathogen, fixed = TRUE)
  status <- ifelse(pathogen == "negative", "virus_negative_symptomatic",
                   ifelse(pathogen == "RSV", "RSV_case", "other_virus_case"))
  include <- !mixed & pathogen != "negative"
  set_rows(
    idx, center = NA_character_, timepoint = "acute_primary", status = status,
    pathogen = pathogen, age_months = numeric_after_colon(raw$characteristics_ch1__repeat_2[idx]),
    sex = normalize_sex(raw$characteristics_ch1__repeat_1[idx]), severity = "hospitalized",
    coinfection = ifelse(mixed, "mixed_virus", "none"),
    batch = paste(after_colon(raw$characteristics_ch1__repeat_4[idx]),
                  after_colon(raw$characteristics_ch1__repeat_5[idx]), sep = ";"),
    tissue = "whole_blood", include = include,
    inclusion_scope = ifelse(include, "exploratory_severity_or_cross_virus", "excluded"),
    exclusion_reason = ifelse(include, NA_character_, ifelse(mixed, "mixed_virus", "virus_negative_symptomatic_hospitalized")),
    needs_manual_review = FALSE,
    metadata_evidence = "sample pathogen+tissue+hospital batch; no healthy controls in this series"
  )

  # GSE103119: only exact RSV single-virus pneumonia and explicit healthy controls.
  idx <- raw$series_accession == "GSE103119"
  viral <- after_colon(raw$characteristics_ch1__repeat_2[idx])
  bacterial <- after_colon(raw$characteristics_ch1__repeat_1[idx])
  condition <- after_colon(raw$characteristics_ch1__repeat_8[idx])
  healthy <- condition == "Healthy Control"
  rsv_single <- condition == "Pneumonia" & viral == "RSV" & bacterial == "NONE"
  include <- healthy | rsv_single
  status <- ifelse(healthy, "healthy_control", ifelse(rsv_single, "RSV_case", "noneligible_pneumonia"))
  subject <- gsub("[[:space:]]+", "_", after_colon(raw$characteristics_ch1[idx]))
  set_rows(
    idx, subject_id = paste0("GSE103119_SUBJ_", subject), center = "Columbus_NCH",
    timepoint = "acute_primary", status = status,
    pathogen = ifelse(healthy, "none", viral),
    age_months = numeric_after_colon(raw$characteristics_ch1__repeat_4[idx]),
    sex = normalize_sex(raw$characteristics_ch1__repeat_5[idx]),
    severity = ifelse(healthy, "healthy", "hospitalized"),
    coinfection = ifelse(rsv_single | healthy, "none", "non_RSV_or_possible_coinfection"),
    tissue = "whole_blood", include = include,
    inclusion_scope = ifelse(include, "exploratory_RSV_single_infection", "excluded"),
    exclusion_reason = ifelse(include, NA_character_, "not_unambiguous_RSV_single_infection_or_healthy"),
    needs_manual_review = FALSE,
    metadata_evidence = "sample name+bacterial organism+viral organism+condition; exact equality required"
  )

  m$meta_eligible_confirmatory <- m$include & m$series_accession %in% c("GSE105450", "GSE103842")
  m$severity_analysis_eligible <- m$include & m$status == "RSV_case" &
    m$severity %in% c("hospitalized", "outpatient") & !m$needs_manual_review
  m <- m[order(m$series_accession, m$sample_id, method = "radix"), , drop = FALSE]
  rownames(m) <- NULL
  m
}

build_flow <- function(m) {
  keys <- ifelse(m$include, paste0("included:", m$inclusion_scope), paste0("excluded:", m$exclusion_reason))
  out <- as.data.frame(table(m$series_accession, keys), stringsAsFactors = FALSE)
  names(out) <- c("series_accession", "decision_reason", "n_samples")
  out <- out[out$n_samples > 0L, ]
  totals <- aggregate(sample_id ~ series_accession, m, length)
  names(totals)[2] <- "official_samples"
  out <- merge(out, totals, by = "series_accession", all.x = TRUE, sort = FALSE)
  out[order(out$series_accession, out$decision_reason, method = "radix"), ]
}

build_missingness <- function(m) {
  fields <- c("subject_id", "center", "platform", "timepoint", "status", "age_months",
              "sex", "severity", "coinfection", "treatment", "batch")
  populations <- list(all_official = rep(TRUE, nrow(m)), included = m$include)
  rows <- list()
  for (series in unique(m$series_accession)) {
    for (population in names(populations)) {
      sel <- m$series_accession == series & populations[[population]]
      for (field in fields) {
        value <- m[[field]][sel]
        missing <- is.na(value) | (is.character(value) & trimws(value) == "")
        rows[[length(rows) + 1L]] <- data.frame(
          series_accession = series, population = population, field = field,
          n_samples = sum(sel), n_missing = sum(missing),
          percent_missing = if (sum(sel)) round(100 * mean(missing), 2) else NA_real_
        )
      }
    }
  }
  do.call(rbind, rows)
}

build_review_packet <- function(raw, manifest) {
  raw <- raw[match(manifest$sample_id, raw$geo_accession), , drop = FALSE]
  stopifnot(identical(manifest$sample_id, raw$geo_accession))
  rule_map <- c(
    GSE105450 = "R01_confirmatory_acute_RSV_or_healthy_deduplicate",
    GSE103842 = "R02_confirmatory_acute_RSV_or_healthy",
    GSE188427 = "R03_substitution_Day1_or_healthy_only",
    GSE38900 = "R04_exploratory_nonNCH_acute_RSV_or_healthy",
    GSE155925 = "R05_exploratory_single_virus_hospitalized",
    GSE103119 = "R06_exploratory_exact_RSV_single_or_healthy"
  )
  subject_evidence_map <- c(
    GSE105450 = "GSM singleton, except exact public-metadata healthy technical pair shares lower-GSM subject key",
    GSE103842 = "one official GSM per child; no repeat relation reported",
    GSE188427 = "title prefix before underscore with -D1/-D30/-D180 removed",
    GSE38900 = "center-namespaced local numeric token from description; RSV visit suffix removed; unknown center falls back to unique GSM",
    GSE155925 = "one official GSM per child; no repeat relation reported",
    GSE103119 = "sample name characteristic"
  )
  char_cols <- c("characteristics_ch1", paste0("characteristics_ch1__repeat_", 1:8))
  char_values <- raw[char_cols]
  names(char_values) <- paste0("raw_characteristic_", seq_along(char_cols))
  packet <- data.frame(
    sample_id = manifest$sample_id,
    series_accession = manifest$series_accession,
    platform = manifest$platform,
    raw_title = raw$title,
    raw_source_name = raw$source_name_ch1,
    raw_description = raw$description,
    char_values,
    rule_id = unname(rule_map[manifest$series_accession]),
    machine_include = manifest$include,
    machine_exclusion_reason = manifest$exclusion_reason,
    machine_status = manifest$status,
    machine_timepoint = manifest$timepoint,
    machine_tissue = manifest$tissue,
    machine_virus = manifest$pathogen,
    machine_severity = manifest$severity,
    metadata_discrepancy = manifest$metadata_discrepancy,
    subject_id = manifest$subject_id,
    subject_id_evidence = unname(subject_evidence_map[manifest$series_accession]),
    review_status = "pending_independent_review",
    reviewer = NA_character_,
    reviewed_at = NA_character_,
    review_decision = NA_character_,
    final_include = NA,
    final_exclusion_reason = NA_character_,
    review_notes = NA_character_,
    packet_schema_version = "task4_review_v1",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  packet
}

review_columns <- c(
  "review_status", "reviewer", "reviewed_at", "review_decision",
  "final_include", "final_exclusion_reason", "review_notes"
)

preserve_existing_reviews <- function(template, path) {
  if (!file.exists(path)) return(template)
  existing <- read.delim(path, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8", na.strings = "")
  if (!identical(names(existing), names(template)) || nrow(existing) != nrow(template) ||
      anyDuplicated(existing$sample_id) || !setequal(existing$sample_id, template$sample_id)) {
    stop("Existing review packet schema/sample IDs differ from the machine template; refusing overwrite.")
  }
  character_review_cols <- setdiff(review_columns, "final_include")
  existing[character_review_cols] <- lapply(existing[character_review_cols], as.character)
  existing$final_include <- as.logical(existing$final_include)
  existing <- existing[match(template$sample_id, existing$sample_id), , drop = FALSE]
  machine_cols <- setdiff(names(template), review_columns)
  canonical <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "<MISSING>"
    x
  }
  same_machine <- vapply(machine_cols, function(name) {
    identical(canonical(existing[[name]]), canonical(template[[name]]))
  }, logical(1))
  if (!all(same_machine)) {
    unsigned_pending <- all(existing$review_status == "pending_independent_review") &&
      all(is.na(trim_na(existing$reviewer))) && all(is.na(trim_na(existing$reviewed_at))) &&
      all(is.na(trim_na(existing$review_decision))) && all(is.na(existing$final_include)) &&
      all(is.na(trim_na(existing$final_exclusion_reason))) && all(is.na(trim_na(existing$review_notes)))
    if (isTRUE(unsigned_pending)) return(template)
    stop("Machine evidence changed after review packet creation: ", paste(machine_cols[!same_machine], collapse = ", "))
  }
  template[review_columns] <- existing[review_columns]
  template
}

validate_independent_review <- function(packet, analyst_id = "sample_audit_analyst") {
  completed <- packet$review_status == "completed_independent_review"
  reviewer_value <- trim_na(packet$reviewer)
  reviewer_ok <- !is.na(reviewer_value) & reviewer_value != analyst_id
  time_value <- trim_na(packet$reviewed_at)
  time_ok <- !is.na(time_value) & grepl(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(Z|[+-][0-9]{4})$",
    time_value
  )
  decision_ok <- packet$review_decision %in% c("agree", "disagree", "unresolved")
  final_known <- !is.na(packet$final_include)
  excluded_reason_ok <- !final_known | packet$final_include |
    !is.na(trim_na(packet$final_exclusion_reason))
  agree_ok <- packet$review_decision != "agree" |
    (final_known & packet$final_include == packet$machine_include &
       (packet$machine_include |
          trim_na(packet$final_exclusion_reason) == trim_na(packet$machine_exclusion_reason)))
  conservative_ok <- !packet$review_decision %in% c("disagree", "unresolved") |
    (final_known & !packet$final_include & !is.na(trim_na(packet$final_exclusion_reason)))
  ok <- completed & reviewer_ok & time_ok & decision_ok & final_known & excluded_reason_ok &
    agree_ok & conservative_ok
  ok[is.na(ok)] <- FALSE
  list(complete = isTRUE(all(ok)), row_ok = ok, n_incomplete = sum(!ok))
}

apply_completed_review <- function(manifest, packet) {
  validation <- validate_independent_review(packet)
  if (!validation$complete) stop("Independent review is incomplete or non-conservative.")
  packet <- packet[match(manifest$sample_id, packet$sample_id), , drop = FALSE]
  manifest$include <- packet$final_include
  manifest$exclusion_reason <- ifelse(packet$final_include, NA_character_, packet$final_exclusion_reason)
  manifest$reviewer <- packet$reviewer
  manifest$reviewed_at <- packet$reviewed_at
  manifest$review_decision <- packet$review_decision
  manifest$review_notes <- packet$review_notes
  manifest$meta_eligible_confirmatory <- manifest$include &
    manifest$series_accession %in% c("GSE105450", "GSE103842")
  manifest$severity_analysis_eligible <- manifest$include & manifest$status == "RSV_case" &
    manifest$severity %in% c("hospitalized", "outpatient") & !manifest$needs_manual_review
  manifest$freeze_status <- "frozen_after_complete_independent_review"
  manifest
}

if (!identical(Sys.getenv("TASK4_SKIP_MAIN"), "1")) {
  raw_path <- file.path("data", "metadata", "all_samples_raw.tsv")
  raw <- read.delim(raw_path, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8")
  manifest <- audit_samples(raw)

  matrix_files <- Sys.glob(file.path("data", "raw", "*", "*series_matrix*.gz"))
  matrix_headers <- lapply(matrix_files, read_series_matrix_header)
  if (!length(matrix_files) || !all(vapply(matrix_headers, `[[`, logical(1), "stopped_before_table"))) {
    stop("Every series matrix must expose a readable metadata header and table boundary.")
  }
  header_audit <- do.call(rbind, Map(function(path, h) data.frame(
    relative_path = path, series_accession = paste(h$accession, collapse = ";"),
    platforms_in_header = paste(h$platforms, collapse = ";"),
    n_sample_ids_in_header = length(h$sample_ids), stopped_before_table = h$stopped_before_table,
    expression_values_read = FALSE
  ), matrix_files, matrix_headers))
  for (series in unique(raw$series_accession)) {
    relevant <- vapply(matrix_headers, function(h) identical(h$accession, series), logical(1))
    header_ids <- unique(unlist(lapply(matrix_headers[relevant], `[[`, "sample_ids"), use.names = FALSE))
    metadata_ids <- raw$geo_accession[raw$series_accession == series]
    if (!setequal(header_ids, metadata_ids)) {
      stop("Series-matrix header sample IDs do not match frozen metadata for ", series)
    }
    quick_path <- file.path("data", "raw", series, paste0(series, "_quick.soft.txt"))
    if (!identical(read_quick_soft_accession(quick_path), series)) {
      stop("Quick SOFT accession mismatch for ", series)
    }
  }

  overlap <- data.frame(
    cohort_a = c("GSE105450", "GSE105450", "GSE38900", "GSE103119"),
    cohort_b = c("GSE188427", "GSE38900", "GSE105450", "GSE105450"),
    overlap_risk = c("high_possible", "center_specific", "center_specific", "possible_NCH_context"),
    evidence = c(
      "same NCH longitudinal study context; participant-level crosswalk unavailable",
      "GSE38900 Columbus is NCH; Dallas/Turku are explicit non-NCH sources",
      "Columbus samples excluded; only Dallas/Turku retained",
      "NCH cohort but exploratory CAP study; no cross-study participant key"
    ),
    frozen_action = c(
      "never include together as independent confirmation; GSE188427 substitutes only",
      "exclude all Columbus samples from GSE38900 exploratory retained set",
      "retained units cannot enter confirmatory meta",
      "exploratory only; never enter confirmatory meta"
    ), stringsAsFactors = FALSE
  )

  manifest$freeze_status <- "preliminary_pending_independent_review"
  review_path <- file.path("data", "clean", "sample_review_packet.tsv")
  review_packet <- preserve_existing_reviews(build_review_packet(raw, manifest), review_path)
  review_validation <- validate_independent_review(review_packet)
  frozen_path <- file.path("data", "clean", "sample_manifest_frozen.tsv")
  if (!review_validation$complete && file.exists(frozen_path)) {
    stop("A stale frozen manifest exists while independent review is incomplete; remove the stale generated artifact before rerun.")
  }

  write_tsv(manifest, file.path("data", "clean", "sample_manifest_preliminary.tsv"))
  write_tsv(review_packet, review_path)

  review_sha <- sha256_file(review_path)
  review_commit <- system2(
    "git", c("log", "-n", "1", "--format=%H", "--", review_path), stdout = TRUE
  )[[1L]]
  working_blob <- system2("git", c("hash-object", review_path), stdout = TRUE)[[1L]]
  committed_blob <- system2(
    "git", c("rev-parse", paste0(review_commit, ":", review_path)), stdout = TRUE
  )[[1L]]
  if (review_validation$complete) {
    stopifnot(
      identical(review_sha, "140239a08915dc88872e7dedfaebe366a00a104ca3c24b936c0861a9229ffdd6"),
      startsWith(review_commit, "fe6dd2a"),
      identical(working_blob, committed_blob),
      all(review_packet$reviewer == "independent_sample_reviewer")
    )
  }

  packet_sha <- data.frame(
    path = review_path, bytes = file.info(review_path)$size,
    sha256 = review_sha, review_status = ifelse(
      review_validation$complete, "complete_independent_review", "pending_independent_review"
    )
  )
  write_tsv(packet_sha, file.path("logs", "checksums", "sample_review_packet_sha256.tsv"))

  if (review_validation$complete) {
    final_manifest <- apply_completed_review(manifest, review_packet)
    write_tsv(final_manifest, frozen_path)
    review_audit <- review_packet[c(
      "sample_id", "series_accession", "platform", "rule_id", "machine_include",
      "machine_exclusion_reason", "machine_status", "machine_severity",
      "metadata_discrepancy", "subject_id", "review_status", "reviewer",
      "reviewed_at", "review_decision", "final_include", "final_exclusion_reason",
      "review_notes", "packet_schema_version"
    )]
    write_tsv(review_audit, file.path("results", "tables", "sample_review_audit.tsv"))
  }

  report_manifest <- if (review_validation$complete) final_manifest else manifest
  write_tsv(build_flow(report_manifest), file.path("results", "tables", "cohort_inclusion_flow.tsv"))
  write_tsv(overlap, file.path("results", "tables", "overlap_risk_register.tsv"))
  write_tsv(build_missingness(report_manifest), file.path("results", "tables", "sample_covariate_missingness.tsv"))
  write_tsv(header_audit, file.path("results", "tables", "series_matrix_header_audit.tsv"))

  confirmatory <- subset(report_manifest, meta_eligible_confirmatory)
  stopifnot(
    !anyDuplicated(report_manifest$sample_id),
    setequal(unique(confirmatory$series_accession), c("GSE105450", "GSE103842")),
    !any(confirmatory$outcome_previously_observed),
    !any(confirmatory$series_accession %in% c("GSE188427", "GSE38900")),
    !any(report_manifest$severity_analysis_eligible & report_manifest$needs_manual_review),
    !any(report_manifest$severity_analysis_eligible & is.na(report_manifest$severity)),
    all(is.na(report_manifest$sex) | report_manifest$sex %in% c("male", "female")),
    nrow(review_packet) == 1031L,
    !anyDuplicated(review_packet$sample_id),
    nrow(report_manifest) == 1031L,
    sum(report_manifest$include) == 611L,
    sum(report_manifest$meta_eligible_confirmatory) == 196L,
    length(unique(confirmatory$subject_id)) == 196L
  )

  log_lines <- c(
    "Task 4: sample inclusion, independent review, and final manifest freeze",
    "Boundary: metadata and series-matrix headers only; expression table rows read=0",
    "Network/downloads=0; source archives deleted=0",
    "Processed-only boundary: no CEL/FASTQ-level QC is claimed",
    sprintf("Official samples=%d; included=%d; excluded=%d", nrow(manifest), sum(manifest$include), sum(!manifest$include)),
    sprintf("Preliminary manifest SHA-256=%s", sha256_file(file.path("data", "clean", "sample_manifest_preliminary.tsv"))),
    sprintf("Independent review packet SHA-256=%s", sha256_file(review_path)),
    sprintf("Independent review packet commit=%s; committed blob matches working packet=%s", review_commit, identical(working_blob, committed_blob)),
    sprintf("Independent review rows pending=%d; final manifest frozen=%s", review_validation$n_incomplete, review_validation$complete),
    if (review_validation$complete) sprintf("Final manifest SHA-256=%s", sha256_file(frozen_path)) else "Final manifest SHA-256=not_created",
    if (review_validation$complete) "Review provenance: reviewer=independent_sample_reviewer; 1031/1031 agree; machine and final inclusion decisions identical" else "Review provenance: pending",
    "GSE105450 technical repeats: five public-metadata pairs; retained assay chosen by lowest GSM accession, independent of expression/outcome",
    "GSE188427 is substitution-only and never simultaneous independent confirmation with GSE105450",
    "GSE188427 metadata conflicts: seven unresolved source_name versus title/characteristics arm records have severity=NA and are ineligible for severity analysis; three Day 1 RSV samples remain in substitution-wide analysis",
    "GSE38900 is previously observed exploratory only; Columbus/NCH samples excluded",
    "GSE38900 subject identity uses center-namespaced local tokens; six cross-center token collisions are separated and twelve Dallas acute/recovery pairs remain linked",
    "No analyst-generated reviewer identity, timestamp, decision, or final inclusion value was written"
  )
  dir.create(file.path("logs", "session_info"), recursive = TRUE, showWarnings = FALSE)
  writeLines(log_lines, file.path("logs", "session_info", "task4_commands.log"), useBytes = TRUE)

  task4_files <- c(
    raw_path, file.path("analysis", "config", "cohorts.yml"),
    file.path("analysis", "R", "03_audit_samples.R"),
    file.path("analysis", "tests", "test_sample_independence.R"),
    file.path("data", "clean", "sample_manifest_preliminary.tsv"),
    review_path,
    file.path("data", "clean", "sample_review_packet_README.md"),
    file.path("results", "tables", "cohort_inclusion_flow.tsv"),
    file.path("results", "tables", "overlap_risk_register.tsv"),
    file.path("results", "tables", "sample_covariate_missingness.tsv"),
    file.path("results", "tables", "series_matrix_header_audit.tsv"),
    file.path("logs", "checksums", "sample_review_packet_sha256.tsv"),
    file.path("logs", "session_info", "task4_commands.log"),
    file.path("logs", "session_info", "task4_gse38900_subject_namespace_anomaly.log"),
    file.path("logs", "session_info", "task4_independent_review.md"),
    file.path("results", "tables", "independent_sample_review_summary.tsv")
  )
  if (review_validation$complete) task4_files <- c(
    task4_files, frozen_path, file.path("results", "tables", "sample_review_audit.tsv")
  )
  checksums <- data.frame(path = task4_files, bytes = file.info(task4_files)$size,
                          sha256 = vapply(task4_files, sha256_file, character(1)))
  write_tsv(checksums, file.path("logs", "checksums", "task4_sha256.tsv"))

  cat(sprintf("Task 4 manifest status=%s; official=%d; included=%d; confirmatory=%d; review pending=%d\n",
              ifelse(review_validation$complete, "frozen", "preliminary"),
              nrow(report_manifest), sum(report_manifest$include),
              sum(report_manifest$meta_eligible_confirmatory), review_validation$n_incomplete))
}
