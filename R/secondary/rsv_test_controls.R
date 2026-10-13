

suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
  library(survival)
  library(broom)
  library(lubridate)
})

SOURCE_MAIN_ANALYSIS <- TRUE

MAIN_ANALYSIS_SCRIPT <- "R/primary/cox_primary.R"

RSV_LAB_FILE <- "kpsc_infant_lb_outcomes_06162026.sas7bdat"

CROSSWALK_FILE <- "linkage_v3_outputs/PROPOSED_crosspull_infant_crosswalk.csv"
if (!file.exists(CROSSWALK_FILE) &&
    file.exists("PROPOSED_crosspull_infant_crosswalk.csv")) {
  CROSSWALK_FILE <- "PROPOSED_crosspull_infant_crosswalk.csv"
}

OUTPUT_DIR <- "rsv_lab_controls_all_abx_pm7day_v10"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

EVENT_RELATIVE_TO_TEST_LOWER <- -7L
EVENT_RELATIVE_TO_TEST_UPPER <- 7L

POSITIVE_EPISODE_GAP_DAYS <- 30L
NEG_AFTER_POS_EXCLUSION_DAYS <- 30L

RESTRICT_TO_RSV_SEASON <- TRUE

USE_3MONTH_CALENDAR_STRATA <- NA

if (SOURCE_MAIN_ANALYSIS) {
  if (!file.exists(MAIN_ANALYSIS_SCRIPT)) {
    stop(
      "Cannot find MAIN_ANALYSIS_SCRIPT: ", MAIN_ANALYSIS_SCRIPT,
      "\nPlace this script in the same analysis directory or edit the path."
    )
  }
  source(MAIN_ANALYSIS_SCRIPT)
}

required_main_objects <- c(
  "cohort_primary", "analytic_events", "blocks_primary", "build_cp",
  "clean_and_restrict", "COVARIATE_TERMS", "FOLLOWUP_MAX"
)
missing_main <- required_main_objects[
  !vapply(required_main_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_main) > 0) {
  stop(
    "Missing canonical antibiotic-analysis objects: ",
    paste(missing_main, collapse = ", "),
    "\nRun/source the complete main antibiotic analysis first."
  )
}

stopifnot(
  n_distinct(cohort_primary$infant_study_id) == nrow(cohort_primary),
  "mom_study_id" %in% names(cohort_primary),
  "index_date" %in% names(cohort_primary)
)

SELECTED_ARI_FLAGS <- c(
  "sup_om_yn",
  "URI_yn",
  "pneumonia_yn",
  "bronchitis_yn"
)
LRTI_FLAGS <- c("pneumonia_yn", "bronchitis_yn")

SELECTED_DX_LABELS <- c(
  sup_om_yn    = "Acute/suppurative otitis media",
  URI_yn       = "Viral upper respiratory infection",
  pneumonia_yn = "Pneumonia",
  bronchitis_yn = "Bronchitis/bronchiolitis"
)

missing_dx_flags <- setdiff(SELECTED_ARI_FLAGS, names(analytic_events))
if (length(missing_dx_flags) > 0) {
  stop("Selected diagnosis flags missing from analytic_events: ",
       paste(missing_dx_flags, collapse = ", "))
}

if (!file.exists(RSV_LAB_FILE)) {
  stop("Cannot find RSV-effectiveness lab file: ", RSV_LAB_FILE)
}

rsv_lab <- read_sas(RSV_LAB_FILE)

required_lab_cols <- c(
  "infant_study_id", "index_date", "outcome_rsvlab_days", "lab_type", "result"
)
missing_lab_cols <- setdiff(required_lab_cols, names(rsv_lab))
if (length(missing_lab_cols) > 0) {
  stop("Required columns missing from RSV lab file: ",
       paste(missing_lab_cols, collapse = ", "))
}

if (!file.exists(CROSSWALK_FILE)) {
  stop("Cannot find cross-pull infant crosswalk: ", CROSSWALK_FILE)
}

crosswalk <- readr::read_csv(CROSSWALK_FILE, show_col_types = FALSE)

required_crosswalk_cols <- c(
  "abx_infant_study_id", "rsv_infant_study_id",
  "index_date_abx", "index_date_rsv",
  "agreement_score", "information_weight"
)
missing_crosswalk_cols <- setdiff(required_crosswalk_cols, names(crosswalk))
if (length(missing_crosswalk_cols) > 0) {
  stop("Crosswalk missing required columns: ",
       paste(missing_crosswalk_cols, collapse = ", "))
}

stopifnot(
  n_distinct(crosswalk$abx_infant_study_id) == nrow(crosswalk),
  n_distinct(crosswalk$rsv_infant_study_id) == nrow(crosswalk)
)

abx_ids <- cohort_primary %>%
  distinct(infant_study_id)

primary_crosswalk <- crosswalk %>%
  semi_join(
    abx_ids,
    by = c("abx_infant_study_id" = "infant_study_id")
  )

primary_linkage_coverage <- tibble(
  metric = c(
    "Infants in canonical antibiotic cohort",
    "Canonical antibiotic-cohort infants with high-confidence cross-pull match",
    "Canonical antibiotic-cohort infants without high-confidence cross-pull match",
    "Percent of canonical antibiotic cohort high-confidence matched"
  ),
  value = c(
    nrow(abx_ids),
    nrow(primary_crosswalk),
    nrow(abx_ids) - nrow(primary_crosswalk),
    100 * nrow(primary_crosswalk) / nrow(abx_ids)
  )
)
print(primary_linkage_coverage, n = Inf)
write_csv(
  primary_linkage_coverage,
  file.path(OUTPUT_DIR, "qa_primary_cohort_crosswalk_coverage.csv")
)

crosswalk_quality_primary <- primary_crosswalk %>%
  summarise(
    n_matches = n(),
    min_agreement_score = min(agreement_score, na.rm = TRUE),
    median_agreement_score = median(agreement_score, na.rm = TRUE),
    max_agreement_score = max(agreement_score, na.rm = TRUE),
    min_information_weight = min(information_weight, na.rm = TRUE),
    median_information_weight = median(information_weight, na.rm = TRUE),
    max_information_weight = max(information_weight, na.rm = TRUE),
    n_date_diff_le_4 = sum(
      abs(as.numeric(as.Date(index_date_rsv) - as.Date(index_date_abx))) <= 4,
      na.rm = TRUE
    )
  )

write_csv(
  crosswalk_quality_primary,
  file.path(OUTPUT_DIR, "qa_primary_cohort_crosswalk_quality.csv")
)

unmatched_primary_ids <- abx_ids %>%
  anti_join(
    primary_crosswalk,
    by = c("infant_study_id" = "abx_infant_study_id")
  )

write_csv(
  unmatched_primary_ids,
  file.path(OUTPUT_DIR, "qa_primary_cohort_unmatched_infants.csv")
)

rsv_lab_crosswalked <- rsv_lab %>%
  inner_join(
    primary_crosswalk %>%
      select(
        abx_infant_study_id,
        rsv_infant_study_id,
        index_date_abx_crosswalk = index_date_abx,
        index_date_rsv_crosswalk = index_date_rsv,
        agreement_score,
        information_weight
      ),
    by = c("infant_study_id" = "rsv_infant_study_id")
  ) %>%
  rename(
    rsv_infant_study_id = infant_study_id,
    infant_study_id = abx_infant_study_id
  )

rsv_lab_ids <- rsv_lab %>%
  distinct(rsv_infant_study_id = infant_study_id)

rsv_lab_linked_ids <- rsv_lab_crosswalked %>%
  distinct(infant_study_id, rsv_infant_study_id)

lab_linkage_summary <- tibble(
  metric = c(
    "Distinct infants in June RSV lab file",
    "Distinct June lab infants linked to canonical antibiotic cohort",
    "Canonical antibiotic-cohort infants with >=1 linked RSV lab record",
    "Canonical antibiotic-cohort infants with high-confidence crosswalk but no RSV lab record",
    "Percent of canonical antibiotic cohort with >=1 linked RSV lab record"
  ),
  value = c(
    nrow(rsv_lab_ids),
    n_distinct(rsv_lab_linked_ids$rsv_infant_study_id),
    n_distinct(rsv_lab_linked_ids$infant_study_id),
    nrow(primary_crosswalk) - n_distinct(rsv_lab_linked_ids$infant_study_id),
    100 * n_distinct(rsv_lab_linked_ids$infant_study_id) / nrow(abx_ids)
  )
)
print(lab_linkage_summary, n = Inf)
write_csv(
  lab_linkage_summary,
  file.path(OUTPUT_DIR, "qa_rsv_lab_crosswalk_linkage_summary.csv")
)

stopifnot(
  rsv_lab_linked_ids %>%
    count(rsv_infant_study_id) %>%
    summarise(ok = all(n == 1)) %>%
    pull(ok)
)

rsv_lab_primary <- rsv_lab_crosswalked
print(rsv_lab_primary %>% count(result, sort = TRUE), n = Inf)
print(rsv_lab_primary %>% count(lab_type, sort = TRUE), n = Inf)

testing <- rsv_lab_primary %>%
  mutate(
    lab_type_clean = case_when(
      grepl("Antigen",  lab_type) ~ "antigen",
      grepl("Subtype A", lab_type) ~ "pcr_a",
      grepl("Subtype B", lab_type) ~ "pcr_b",
      TRUE                         ~ "other"
    ),
    detected = result %in% c("Positive", "Detected")
  ) %>%
  group_by(infant_study_id, index_date, outcome_rsvlab_days) %>%
  summarise(
    rsv_pos   = any(detected),
    rsv_a_pos = any(detected & lab_type_clean == "pcr_a"),
    rsv_b_pos = any(detected & lab_type_clean == "pcr_b"),
    rsv_neg   = all(!detected),
    n_lab_rows = n(),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(outcome_rsvlab_days),
    outcome_rsvlab_days > 0,
    outcome_rsvlab_days <= FOLLOWUP_MAX
  )

stopifnot(!any(testing$rsv_pos & testing$rsv_neg))

lab_collapse_summary <- tibble(
  metric = c(
    "Raw RSV lab rows in linked primary cohort",
    "Distinct RSV test-days, days 1-followup",
    "RSV-positive test-days",
    "RSV-negative test-days"
  ),
  n = c(
    nrow(rsv_lab_primary),
    nrow(testing),
    sum(testing$rsv_pos),
    sum(testing$rsv_neg)
  )
)
print(lab_collapse_summary, n = Inf)
write_csv(lab_collapse_summary, file.path(OUTPUT_DIR, "qa_lab_collapse_summary.csv"))

cases_raw <- testing %>%
  filter(rsv_pos) %>%
  arrange(infant_study_id, outcome_rsvlab_days)

cases <- cases_raw %>%
  group_by(infant_study_id) %>%
  filter({
    keep <- logical(n())
    last_kept <- -Inf
    for (i in seq_len(n())) {
      if (outcome_rsvlab_days[i] - last_kept > POSITIVE_EPISODE_GAP_DAYS) {
        keep[i] <- TRUE
        last_kept <- outcome_rsvlab_days[i]
      }
    }
    keep
  }) %>%
  ungroup()

controls <- testing %>%
  filter(rsv_neg) %>%
  left_join(
    cases %>% select(infant_study_id, pos_days = outcome_rsvlab_days),
    by = "infant_study_id",
    relationship = "many-to-many"
  ) %>%
  group_by(infant_study_id, outcome_rsvlab_days) %>%
  filter(!any(
    !is.na(pos_days) &
      outcome_rsvlab_days >= pos_days &
      outcome_rsvlab_days <= pos_days + NEG_AFTER_POS_EXCLUSION_DAYS
  )) %>%
  ungroup() %>%
  select(-pos_days) %>%
  distinct()

test_episode_reconciliation <- tibble(
  stage = c(
    "Positive test-days before 30-day episode collapse",
    "Positive episodes retained",
    "Positive test-days dropped as repeat positives",
    "Negative test-days before exclusions",
    "Negative test-days after 0-30d post-positive exclusion"
  ),
  n = c(
    nrow(cases_raw),
    nrow(cases),
    nrow(cases_raw) - nrow(cases),
    sum(testing$rsv_neg),
    nrow(controls)
  )
)
print(test_episode_reconciliation, n = Inf)
write_csv(test_episode_reconciliation,
          file.path(OUTPUT_DIR, "qa_test_episode_reconciliation.csv"))

clean_test_events <- bind_rows(
  cases %>% transmute(
    infant_study_id,
    test_day = outcome_rsvlab_days,
    lab_status = "RSV-positive"
  ),
  controls %>% transmute(
    infant_study_id,
    test_day = outcome_rsvlab_days,
    lab_status = "RSV-negative"
  )
) %>%
  distinct(infant_study_id, test_day, lab_status)

stopifnot(
  !anyDuplicated(clean_test_events[c("infant_study_id", "test_day", "lab_status")])
)

required_event_cols <- c(
  "infant_study_id", "abx_days", "analyzable_event_primary", SELECTED_ARI_FLAGS
)
stopifnot(all(required_event_cols %in% names(analytic_events)))

primary_antibiotic_records <- analytic_events %>%
  filter(
    analyzable_event_primary == TRUE,
    !is.na(abx_days),
    abx_days > 0,
    abx_days <= FOLLOWUP_MAX
  ) %>%
  semi_join(abx_ids, by = "infant_study_id")

abx_day_dx <- primary_antibiotic_records %>%
  group_by(infant_study_id, abx_days) %>%
  summarise(
    across(all_of(SELECTED_ARI_FLAGS), ~ as.integer(any(.x == 1, na.rm = TRUE))),
    .groups = "drop"
  ) %>%
  mutate(
    selected_ari = as.integer(if_any(all_of(SELECTED_ARI_FLAGS), ~ .x == 1)),
    lrti = as.integer(if_any(all_of(LRTI_FLAGS), ~ .x == 1))
  )

all_abx_days <- primary_antibiotic_records %>%
  distinct(infant_study_id, abx_days)

stopifnot(
  nrow(abx_day_dx) == nrow(all_abx_days),
  all(abx_day_dx$lrti <= abx_day_dx$selected_ari)
)

individual_selected_union <- bind_rows(
  lapply(SELECTED_ARI_FLAGS, function(flag) {
    abx_day_dx %>%
      filter(.data[[flag]] == 1) %>%
      select(infant_study_id, abx_days)
  })
) %>%
  distinct()

selected_ari_days <- abx_day_dx %>%
  filter(selected_ari == 1) %>%
  select(infant_study_id, abx_days)

stopifnot(
  nrow(anti_join(selected_ari_days, individual_selected_union,
                 by = c("infant_study_id", "abx_days"))) == 0,
  nrow(anti_join(individual_selected_union, selected_ari_days,
                 by = c("infant_study_id", "abx_days"))) == 0
)

abx_dx_reconciliation <- tibble(
  metric = c(
    "Systemic antibiotic infant-days, days 1-followup",
    "Selected ARI-associated antibiotic infant-days",
    "LRTI-associated antibiotic infant-days",
    "Acute/suppurative otitis media-associated antibiotic infant-days",
    "Viral URI-associated antibiotic infant-days",
    "Pneumonia-associated antibiotic infant-days",
    "Bronchitis/bronchiolitis-associated antibiotic infant-days"
  ),
  n = c(
    nrow(all_abx_days),
    sum(abx_day_dx$selected_ari),
    sum(abx_day_dx$lrti),
    sum(abx_day_dx$sup_om_yn),
    sum(abx_day_dx$URI_yn),
    sum(abx_day_dx$pneumonia_yn),
    sum(abx_day_dx$bronchitis_yn)
  )
)
print(abx_dx_reconciliation, n = Inf)
write_csv(abx_dx_reconciliation,
          file.path(OUTPUT_DIR, "qa_antibiotic_dx_event_reconciliation.csv"))

linked_event_tests <- abx_day_dx %>%
  inner_join(clean_test_events, by = "infant_study_id", relationship = "many-to-many") %>%
  mutate(event_relative_to_test = abx_days - test_day) %>%
  filter(
    event_relative_to_test >= EVENT_RELATIVE_TO_TEST_LOWER,
    event_relative_to_test <= EVENT_RELATIVE_TO_TEST_UPPER
  )

status_per_event <- linked_event_tests %>%
  group_by(infant_study_id, abx_days) %>%
  summarise(
    any_positive_test = any(lab_status == "RSV-positive"),
    any_negative_test = any(lab_status == "RSV-negative"),
    n_linked_positive_tests = sum(lab_status == "RSV-positive"),
    n_linked_negative_tests = sum(lab_status == "RSV-negative"),
    .groups = "drop"
  ) %>%
  mutate(
    rsv_event_status = case_when(
      any_positive_test ~ "RSV-positive",
      !any_positive_test & any_negative_test ~ "RSV-negative",
      TRUE ~ NA_character_
    )
  )

abx_dx_tested <- abx_day_dx %>%
  left_join(status_per_event, by = c("infant_study_id", "abx_days"))

pos_event_keys <- abx_dx_tested %>%
  filter(rsv_event_status == "RSV-positive") %>%
  distinct(infant_study_id, abx_days)
neg_event_keys <- abx_dx_tested %>%
  filter(rsv_event_status == "RSV-negative") %>%
  distinct(infant_study_id, abx_days)

stopifnot(
  nrow(inner_join(pos_event_keys, neg_event_keys,
                  by = c("infant_study_id", "abx_days"))) == 0
)

linkage_reconciliation <- tibble(
  metric = c(
    "All systemic antibiotic infant-days",
    "All systemic antibiotic infant-days with any linked clean RSV test",
    "RSV-positive antibiotic infant-days",
    "RSV-negative antibiotic infant-days",
    "Events with both positive and negative linked tests (classified positive)",
    "All systemic antibiotic infant-days without linked clean RSV test"
  ),
  n = c(
    nrow(abx_day_dx),
    nrow(status_per_event),
    nrow(pos_event_keys),
    nrow(neg_event_keys),
    sum(status_per_event$any_positive_test & status_per_event$any_negative_test),
    nrow(abx_day_dx) - nrow(status_per_event)
  )
)
print(linkage_reconciliation, n = Inf)
write_csv(linkage_reconciliation,
          file.path(OUTPUT_DIR, "qa_event_test_linkage_reconciliation.csv"))

linkage_distance_summary <- linked_event_tests %>%
  count(lab_status, event_relative_to_test, name = "n_event_test_links") %>%
  arrange(lab_status, event_relative_to_test)
write_csv(linkage_distance_summary,
          file.path(OUTPUT_DIR, "qa_event_test_linkage_distance.csv"))

linkage_window_definition <- tibble(
  parameter = c(
    "event_relative_to_test_lower_days",
    "event_relative_to_test_upper_days",
    "window_width_days_inclusive",
    "positive_status_precedence"
  ),
  value = c(
    as.character(EVENT_RELATIVE_TO_TEST_LOWER),
    as.character(EVENT_RELATIVE_TO_TEST_UPPER),
    as.character(EVENT_RELATIVE_TO_TEST_UPPER -
                   EVENT_RELATIVE_TO_TEST_LOWER + 1L),
    "Yes"
  )
)
write_csv(
  linkage_window_definition,
  file.path(OUTPUT_DIR, "qa_event_test_linkage_window_definition.csv")
)

ALL_ABX_CONTROL_EVENT_SETS <- list(
  "RSV-positive antibiotic use" = abx_dx_tested %>%
    filter(rsv_event_status == "RSV-positive") %>%
    distinct(infant_study_id, abx_days),
  "RSV-negative antibiotic use" = abx_dx_tested %>%
    filter(rsv_event_status == "RSV-negative") %>%
    distinct(infant_study_id, abx_days)
)

all_abx_control_event_counts_raw <- imap_dfr(
  ALL_ABX_CONTROL_EVENT_SETS,
  ~ tibble(
    outcome = .y,
    n_events_raw = nrow(.x),
    n_infants_raw = n_distinct(.x$infant_study_id)
  )
)

write_csv(
  all_abx_control_event_counts_raw,
  file.path(OUTPUT_DIR, "qa_all_antibiotic_control_event_sets.csv")
)

make_status_event_set <- function(status, flag = "selected_ari") {
  if (!flag %in% c("selected_ari", "lrti", SELECTED_ARI_FLAGS)) {
    stop("Unknown outcome flag: ", flag)
  }

  abx_dx_tested %>%
    filter(
      rsv_event_status == status,
      .data[[flag]] == 1
    ) %>%
    distinct(infant_study_id, abx_days)
}

CONTROL_EVENT_SETS <- list(
  "RSV-positive selected ARI" = make_status_event_set("RSV-positive", "selected_ari"),
  "RSV-negative selected ARI" = make_status_event_set("RSV-negative", "selected_ari"),
  "RSV-positive LRTI" = make_status_event_set("RSV-positive", "lrti"),
  "RSV-negative LRTI" = make_status_event_set("RSV-negative", "lrti")
)

for (flag in SELECTED_ARI_FLAGS) {
  label <- SELECTED_DX_LABELS[[flag]]
  CONTROL_EVENT_SETS[[paste0("RSV-positive: ", label)]] <-
    make_status_event_set("RSV-positive", flag)
  CONTROL_EVENT_SETS[[paste0("RSV-negative: ", label)]] <-
    make_status_event_set("RSV-negative", flag)
}

qa_event_funnel <- tibble(
  stage = c(
    "All systemic antibiotic infant-days",
    "All antibiotic infant-days with qualifying RSV test +/-7 days",
    "RSV-positive antibiotic infant-days",
    "RSV-negative antibiotic infant-days",
    "RSV-positive selected-diagnosis antibiotic infant-days",
    "RSV-negative selected-diagnosis antibiotic infant-days"
  ),
  n_events = c(
    nrow(abx_day_dx),
    sum(!is.na(abx_dx_tested$rsv_event_status)),
    nrow(ALL_ABX_CONTROL_EVENT_SETS[["RSV-positive antibiotic use"]]),
    nrow(ALL_ABX_CONTROL_EVENT_SETS[["RSV-negative antibiotic use"]]),
    nrow(CONTROL_EVENT_SETS[["RSV-positive selected ARI"]]),
    nrow(CONTROL_EVENT_SETS[["RSV-negative selected ARI"]])
  ),
  n_infants = c(
    n_distinct(abx_day_dx$infant_study_id),
    n_distinct(abx_dx_tested$infant_study_id[!is.na(abx_dx_tested$rsv_event_status)]),
    n_distinct(ALL_ABX_CONTROL_EVENT_SETS[["RSV-positive antibiotic use"]]$infant_study_id),
    n_distinct(ALL_ABX_CONTROL_EVENT_SETS[["RSV-negative antibiotic use"]]$infant_study_id),
    n_distinct(CONTROL_EVENT_SETS[["RSV-positive selected ARI"]]$infant_study_id),
    n_distinct(CONTROL_EVENT_SETS[["RSV-negative selected ARI"]]$infant_study_id)
  )
)

write_csv(qa_event_funnel, file.path(OUTPUT_DIR, "qa_all_antibiotic_event_funnel.csv"))
print(qa_event_funnel, n = Inf)

n_all_pos <- nrow(ALL_ABX_CONTROL_EVENT_SETS[["RSV-positive antibiotic use"]])
n_all_neg <- nrow(ALL_ABX_CONTROL_EVENT_SETS[["RSV-negative antibiotic use"]])
n_sel_pos <- nrow(CONTROL_EVENT_SETS[["RSV-positive selected ARI"]])
n_sel_neg <- nrow(CONTROL_EVENT_SETS[["RSV-negative selected ARI"]])

stopifnot(
  length(n_all_pos) == 1L, length(n_all_neg) == 1L,
  length(n_sel_pos) == 1L, length(n_sel_neg) == 1L,
  !anyNA(c(n_all_pos, n_all_neg, n_sel_pos, n_sel_neg)),
  n_all_pos >= n_sel_pos, n_all_neg >= n_sel_neg
)

if (isTRUE(n_all_pos == n_sel_pos) || isTRUE(n_all_neg == n_sel_neg)) {
  warning(
    "At least one all-antibiotic control count equals its selected-diagnosis count. ",
    "Review qa_all_antibiotic_event_funnel.csv to determine whether this is genuine."
  )
}

diagnosis_restricted_selected_counts_raw <- tibble(
  outcome = c(
    "RSV-positive selected-diagnosis antibiotic use",
    "RSV-negative selected-diagnosis antibiotic use"
  ),
  n_events_raw = c(n_sel_pos, n_sel_neg),
  n_infants_raw = c(
    n_distinct(CONTROL_EVENT_SETS[["RSV-positive selected ARI"]]$infant_study_id),
    n_distinct(CONTROL_EVENT_SETS[["RSV-negative selected ARI"]]$infant_study_id)
  )
)

write_csv(
  bind_rows(
    all_abx_control_event_counts_raw %>%
      mutate(definition = "All antibiotic use within +/-7 days of qualifying RSV test"),
    diagnosis_restricted_selected_counts_raw %>%
      mutate(definition = "Selected-diagnosis antibiotic use within +/-7 days of qualifying RSV test")
  ),
  file.path(OUTPUT_DIR, "qa_all_vs_diagnosis_restricted_event_counts.csv")
)

event_set_qa <- imap_dfr(CONTROL_EVENT_SETS, function(events_df, outcome_name) {
  duplicates <- nrow(events_df) - nrow(distinct(events_df, infant_study_id, abx_days))
  outside_fu <- sum(events_df$abx_days <= 0 | events_df$abx_days > FOLLOWUP_MAX, na.rm = TRUE)
  not_in_selected_ari <- nrow(anti_join(
    events_df,
    selected_ari_days,
    by = c("infant_study_id", "abx_days")
  ))

  tibble(
    outcome = outcome_name,
    n_events = nrow(events_df),
    n_infants_with_event = n_distinct(events_df$infant_study_id),
    duplicate_infant_days = duplicates,
    events_outside_followup = outside_fu,
    events_not_in_selected_ari = not_in_selected_ari,
    min_event_day = ifelse(nrow(events_df) > 0, min(events_df$abx_days), NA),
    max_event_day = ifelse(nrow(events_df) > 0, max(events_df$abx_days), NA)
  )
})
print(event_set_qa, n = Inf)
write_csv(event_set_qa, file.path(OUTPUT_DIR, "qa_control_event_sets.csv"))

stopifnot(
  all(event_set_qa$duplicate_infant_days == 0),
  all(event_set_qa$events_outside_followup == 0),
  all(event_set_qa$events_not_in_selected_ari == 0)
)

for (status in c("RSV-positive", "RSV-negative")) {
  ari_name <- paste(status, "selected ARI")
  lrti_name <- paste(status, "LRTI")
  stopifnot(
    nrow(anti_join(
      CONTROL_EVENT_SETS[[lrti_name]],
      CONTROL_EVENT_SETS[[ari_name]],
      by = c("infant_study_id", "abx_days")
    )) == 0
  )
}

add_calendar_3month <- function(df) {
  df %>%
    mutate(
      calendar_3month = factor(paste0(
        format(calendar_date_start, "%Y"), "_Q",
        (as.integer(format(calendar_date_start, "%m")) - 1L) %/% 3L + 1L
      ))
    )
}

control_formula_1month <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_month) +",
  COVARIATE_TERMS
))

control_formula_3month <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_3month) +",
  COVARIATE_TERMS
))

model_arm_counts <- function(md, outcome_name) {
  md %>%
    group_by(exposure_arm) %>%
    summarise(
      n_infants_contributing = n_distinct(infant_study_id),
      person_days = sum(tstop - tstart),
      person_years = person_days / 365.25,
      n_events = sum(event),
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_name, .before = 1)
}

format_exposure <- function(x) {
  case_when(
    x == "d_neither" ~ "Neither",
    x == "a_maternal_vax_only" ~ "Maternal RSV vaccine",
    x == "b_nirsevimab_only" ~ "Nirsevimab",
    x == "c_both" ~ "Both",
    TRUE ~ as.character(x)
  )
}

fit_control_outcome <- function(
    events_df,
    outcome_name,
    restrict_to_season = RESTRICT_TO_RSV_SEASON,
    use_3month_strata = FALSE
) {

  events_df <- events_df %>%
    filter(!is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX) %>%
    semi_join(abx_ids, by = "infant_study_id") %>%
    distinct(infant_study_id, abx_days)

  cp <- build_cp(outcome_events = events_df, blocks = blocks_primary)

  md <- clean_and_restrict(
    cp,
    restrict_to_season = restrict_to_season
  )

  if (use_3month_strata) {
    md <- add_calendar_3month(md)
    model_formula <- control_formula_3month
    strata_specification <- "3-month calendar strata"
  } else {
    model_formula <- control_formula_1month
    strata_specification <- "1-month calendar strata"
  }

  if ("c_both" %in% as.character(unique(md$exposure_arm))) {
    n_both_events <- md %>%
      filter(exposure_arm == "c_both") %>%
      summarise(n = sum(event)) %>%
      pull(n)
    if (n_both_events > 0) {
      stop(outcome_name, ": observed events in c_both despite never-both cohort.")
    }
    md <- md %>%
      filter(exposure_arm != "c_both") %>%
      mutate(exposure_arm = droplevels(exposure_arm))
  }

  counts <- model_arm_counts(md, outcome_name)

  model_n_events <- sum(md$event)
  arm_n_events <- sum(counts$n_events)
  model_person_days <- sum(md$tstop - md$tstart)
  arm_person_days <- sum(counts$person_days)

  if (model_n_events != arm_n_events) {
    stop(outcome_name, ": model events do not equal exposure-arm event sum.")
  }
  if (abs(model_person_days - arm_person_days) > 1e-8) {
    stop(outcome_name, ": model person-days do not equal exposure-arm sum.")
  }

  fit <- tryCatch(
    coxph(
      model_formula,
      data = md,
      cluster = mom_study_id,
      ties = "efron",
      model = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      model_data = md,
      counts = counts,
      results = tibble(
        outcome = outcome_name,
        strata_specification = strata_specification,
        exposure = NA_character_,
        HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
        p_value = NA_real_, fit_status = paste("error:", conditionMessage(fit)),
        n_events_raw = nrow(events_df),
        n_events_model = model_n_events
      ),
      reconciliation = tibble(
        outcome = outcome_name,
        strata_specification = strata_specification,
        n_events_raw = nrow(events_df),
        n_events_model = model_n_events,
        n_events_from_arms = arm_n_events,
        person_days_model = model_person_days,
        person_days_from_arms = arm_person_days
      )
    ))
  }

  td <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(str_starts(term, "exposure_arm")) %>%
    transmute(
      outcome = outcome_name,
      strata_specification = strata_specification,
      exposure_term = term,
      exposure_arm = str_remove(term, "^exposure_arm"),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      VE_pct = 100 * (1 - estimate),
      VE_CI_low = 100 * (1 - conf.high),
      VE_CI_high = 100 * (1 - conf.low),
      p_value = p.value,
      fit_status = "fit",
      n_events_raw = nrow(events_df),
      n_events_model = model_n_events
    ) %>%
    left_join(
      counts %>%
        select(exposure_arm, n_infants_contributing, person_days,
               person_years, n_events),
      by = "exposure_arm"
    ) %>%
    mutate(exposure = format_exposure(exposure_arm))

  ref <- counts %>%
    filter(exposure_arm == "d_neither") %>%
    transmute(
      outcome = outcome_name,
      strata_specification = strata_specification,
      exposure_term = "Reference",
      exposure_arm,
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
      p_value = NA_real_, fit_status = "fit",
      n_events_raw = nrow(events_df),
      n_events_model = model_n_events,
      n_infants_contributing,
      person_days,
      person_years,
      n_events,
      exposure = "Neither"
    )

  results <- bind_rows(td, ref) %>%
    mutate(
      exposure = factor(
        exposure,
        levels = c("Neither", "Maternal RSV vaccine", "Nirsevimab")
      )
    ) %>%
    arrange(exposure)

  list(
    model_data = md,
    counts = counts,
    results = results,
    reconciliation = tibble(
      outcome = outcome_name,
      strata_specification = strata_specification,
      n_events_raw = nrow(events_df),
      n_events_model = model_n_events,
      n_events_from_arms = arm_n_events,
      person_days_model = model_person_days,
      person_days_from_arms = arm_person_days
    )
  )
}

primary_neg_season <- fit_control_outcome(
  ALL_ABX_CONTROL_EVENT_SETS[["RSV-negative antibiotic use"]],
  "RSV-negative test-associated antibiotic use | RSV season",
  restrict_to_season = TRUE,
  use_3month_strata = FALSE
)

positive_control_season <- fit_control_outcome(
  ALL_ABX_CONTROL_EVENT_SETS[["RSV-positive antibiotic use"]],
  "RSV-positive test-associated antibiotic use | RSV season",
  restrict_to_season = TRUE,
  use_3month_strata = FALSE
)

negative_full_year <- fit_control_outcome(
  ALL_ABX_CONTROL_EVENT_SETS[["RSV-negative antibiotic use"]],
  "RSV-negative test-associated antibiotic use | full year",
  restrict_to_season = FALSE,
  use_3month_strata = FALSE
)

positive_full_year <- fit_control_outcome(
  ALL_ABX_CONTROL_EVENT_SETS[["RSV-positive antibiotic use"]],
  "RSV-positive test-associated antibiotic use | full year",
  restrict_to_season = FALSE,
  use_3month_strata = FALSE
)

PRESPECIFIED_MODELS <- list(
  "1_primary_negative_control_rsv_season" = primary_neg_season,
  "2_positive_control_rsv_season" = positive_control_season,
  "3a_full_year_negative_sensitivity" = negative_full_year,
  "3b_full_year_positive_sensitivity" = positive_full_year
)

prespecified_results <- imap_dfr(
  PRESPECIFIED_MODELS,
  ~ .x$results %>% mutate(analysis = .y, .before = 1)
)
prespecified_counts <- imap_dfr(
  PRESPECIFIED_MODELS,
  ~ .x$counts %>% mutate(analysis = .y, .before = 1)
)
prespecified_reconciliation <- imap_dfr(
  PRESPECIFIED_MODELS,
  ~ .x$reconciliation %>% mutate(analysis = .y, .before = 1)
)

write_csv(prespecified_results,
          file.path(OUTPUT_DIR, "prespecified_all_abx_results.csv"))
write_csv(prespecified_counts,
          file.path(OUTPUT_DIR, "prespecified_all_abx_arm_counts.csv"))
write_csv(prespecified_reconciliation,
          file.path(OUTPUT_DIR, "prespecified_all_abx_model_reconciliation.csv"))

prespecified_manuscript_table <- prespecified_results %>%
  transmute(
    analysis, outcome, strata_specification,
    exposure = as.character(exposure),
    n_events, n_infants_contributing, person_years,
    HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high,
    p_value, fit_status
  )

write_csv(
  prespecified_manuscript_table,
  file.path(OUTPUT_DIR, "prespecified_all_abx_manuscript_table.csv")
)
print(prespecified_results, n = Inf)

model_objects <- imap(
  CONTROL_EVENT_SETS,
  ~ fit_control_outcome(
      .x,
      .y,
      restrict_to_season = TRUE,
      use_3month_strata = TRUE
    )
)

control_results <- map_dfr(model_objects, "results")
control_counts <- map_dfr(model_objects, "counts")
model_reconciliation <- map_dfr(model_objects, "reconciliation")

prespecified_model_qa <- prespecified_reconciliation %>%
  mutate(
    event_counts_reconcile =
      n_events_model == n_events_from_arms,
    person_time_reconciles =
      abs(person_days_model - person_days_from_arms) < 1e-8
  )

write_csv(
  prespecified_model_qa,
  file.path(OUTPUT_DIR, "qa_prespecified_all_abx_models.csv")
)

if (!all(prespecified_model_qa$event_counts_reconcile &
         prespecified_model_qa$person_time_reconciles)) {
  stop("One or more prespecified all-antibiotic control models failed count/person-time reconciliation.")
}

model_reconciliation <- model_reconciliation %>%
  mutate(
    event_difference = n_events_model - n_events_from_arms,
    person_day_difference = person_days_model - person_days_from_arms,
    pct_raw_events_entering_model = if_else(
      n_events_raw > 0,
      100 * n_events_model / n_events_raw,
      NA_real_
    )
  )
print(model_reconciliation, n = Inf)
write_csv(model_reconciliation,
          file.path(OUTPUT_DIR, "qa_model_reconciliation.csv"))

stopifnot(
  all(model_reconciliation$event_difference == 0),
  all(abs(model_reconciliation$person_day_difference) < 1e-8)
)

cp_integrity <- imap_dfr(model_objects, function(obj, outcome_name) {
  md <- obj$model_data
  tibble(
    outcome = outcome_name,
    n_rows = nrow(md),
    duplicate_interval_rows = nrow(md) - nrow(distinct(
      md, infant_study_id, tstart, tstop, exposure_arm
    )),
    n_infants = n_distinct(md$infant_study_id),
    n_events = sum(md$event),
    min_tstart = ifelse(nrow(md) > 0, min(md$tstart), NA),
    max_tstop = ifelse(nrow(md) > 0, max(md$tstop), NA)
  )
})
print(cp_integrity, n = Inf)
write_csv(cp_integrity, file.path(OUTPUT_DIR, "qa_counting_process_integrity.csv"))

stopifnot(
  all(cp_integrity$duplicate_interval_rows == 0),
  all(cp_integrity$min_tstart >= 0, na.rm = TRUE),
  all(cp_integrity$max_tstop <= FOLLOWUP_MAX, na.rm = TRUE)
)

control_results <- control_results %>%
  mutate(
    control_role = case_when(
      str_starts(outcome, "RSV-positive") ~ "Positive control",
      str_starts(outcome, "RSV-negative") ~ "Negative control",
      TRUE ~ NA_character_
    )
  )

primary_outcome_names <- c(
  "RSV-positive selected ARI",
  "RSV-negative selected ARI",
  "RSV-positive LRTI",
  "RSV-negative LRTI"
)

primary_control_results <- control_results %>%
  filter(outcome %in% primary_outcome_names) %>%
  mutate(
    outcome = factor(outcome, levels = primary_outcome_names)
  ) %>%
  arrange(outcome, exposure)
print(
  primary_control_results %>%
    select(
      outcome, control_role, exposure,
      n_infants_contributing, person_years, n_events,
      HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high,
      p_value, fit_status
    ),
  n = Inf,
  width = Inf
)

write_csv(control_results,
          file.path(OUTPUT_DIR, "rsv_positive_negative_control_results_all.csv"))
write_csv(primary_control_results,
          file.path(OUTPUT_DIR, "rsv_positive_negative_control_results_primary.csv"))
write_csv(control_counts,
          file.path(OUTPUT_DIR, "rsv_positive_negative_control_arm_counts.csv"))

manuscript_table <- primary_control_results %>%
  mutate(
    `HR (95% CI)` = if_else(
      exposure == "Neither",
      "Reference",
      sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high)
    ),
    `Effectiveness, % (95% CI)` = if_else(
      exposure == "Neither",
      "Reference",
      sprintf("%.1f (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high)
    )
  ) %>%
  select(
    Outcome = outcome,
    `Control role` = control_role,
    Exposure = exposure,
    `Infants contributing person-time` = n_infants_contributing,
    `Person-years` = person_years,
    Events = n_events,
    `HR (95% CI)`,
    `Effectiveness, % (95% CI)`
  )

write_csv(manuscript_table,
          file.path(OUTPUT_DIR, "rsv_positive_negative_control_manuscript_table.csv"))

final_qa <- tibble(
  check = c(
    "Primary cohort infant IDs are unique",
    "High-confidence crosswalk index-date differences are within +/-4 days",
    "Collapsed test-days cannot be simultaneously positive and negative",
    "Positive and negative antibiotic event sets are mutually exclusive",
    "LRTI event sets are subsets of selected ARI event sets",
    "All event sets contain unique infant-days",
    "All event sets are within follow-up",
    "All model event totals reconcile to exposure-arm totals",
    "All model person-time totals reconcile to exposure-arm totals",
    "Counting-process interval keys are unique"
  ),
  passed = c(
    n_distinct(cohort_primary$infant_study_id) == nrow(cohort_primary),
    nrow(primary_crosswalk) > 0 &&
      all(
        abs(
          as.numeric(
            as.Date(primary_crosswalk$index_date_rsv) -
              as.Date(primary_crosswalk$index_date_abx)
          )
        ) <= 4,
        na.rm = TRUE
      ),
    !any(testing$rsv_pos & testing$rsv_neg),
    nrow(inner_join(pos_event_keys, neg_event_keys,
                    by = c("infant_study_id", "abx_days"))) == 0,
    TRUE,
    all(event_set_qa$duplicate_infant_days == 0),
    all(event_set_qa$events_outside_followup == 0),
    all(model_reconciliation$event_difference == 0),
    all(abs(model_reconciliation$person_day_difference) < 1e-8),
    all(cp_integrity$duplicate_interval_rows == 0)
  )
)
print(final_qa, n = Inf)
write_csv(final_qa, file.path(OUTPUT_DIR, "qa_final_summary.csv"))

if (!all(final_qa$passed)) {
  stop("One or more final QA checks failed. Review the QA CSV files before reporting results.")
}
