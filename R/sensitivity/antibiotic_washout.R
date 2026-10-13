

library(tidyverse)
library(survival)
library(broom)

stopifnot(
  exists("events_primary"),
  exists("apply_washout"),
  exists("build_cp"),
  exists("clean_and_restrict"),
  exists("fit_ag"),
  exists("tidy_ve"),
  exists("cohort_primary"),
  exists("arm_counts"),
  exists("check_arm_counts"),
  exists("FOLLOWUP_MAX")
)

WASHOUT_SENS_DAYS <- c(14L, 30L)

base_events <- events_primary$all_cause %>%
  distinct(infant_study_id, abx_days) %>%
  arrange(infant_study_id, abx_days)

stopifnot(
  !anyDuplicated(base_events[c("infant_study_id", "abx_days")]),
  all(base_events$abx_days > 0, na.rm = TRUE)
)

washout_events <- setNames(
  lapply(WASHOUT_SENS_DAYS, function(w) apply_washout(base_events, w)),
  paste0("washout_", WASHOUT_SENS_DAYS, "d")
)

event_definition_counts <- tibble(
  definition = c("Primary: unique antibiotic-use days",
                 paste0(WASHOUT_SENS_DAYS, "-day washout")),
  washout_days = c(0L, WASHOUT_SENS_DAYS),
  n_events = c(nrow(base_events), map_int(washout_events, nrow)),
  n_infants_with_event = c(
    n_distinct(base_events$infant_study_id),
    map_int(washout_events, ~ n_distinct(.x$infant_study_id))
  )
)

stopifnot(
  event_definition_counts$n_events[1] >= event_definition_counts$n_events[2],
  event_definition_counts$n_events[2] >= event_definition_counts$n_events[3]
)
print(event_definition_counts, n = Inf)

gap_distribution <- base_events %>%
  group_by(infant_study_id) %>%
  arrange(abx_days, .by_group = TRUE) %>%
  mutate(gap_from_previous_abx_day = abx_days - lag(abx_days)) %>%
  ungroup() %>%
  summarise(
    adjacent_pairs = sum(!is.na(gap_from_previous_abx_day)),
    gap_le_14 = sum(gap_from_previous_abx_day <= 14, na.rm = TRUE),
    gap_15_30 = sum(gap_from_previous_abx_day > 14 &
                      gap_from_previous_abx_day <= 30, na.rm = TRUE),
    gap_gt_30 = sum(gap_from_previous_abx_day > 30, na.rm = TRUE)
  )
print(gap_distribution)

cp_washout <- map(washout_events, build_cp)

make_period_datasets <- function(cp) {
  list(
    rsv_season = clean_and_restrict(
      cp,
      restrict_to_season = TRUE,
      keep_out_of_season_only = FALSE
    ),
    full_year = clean_and_restrict(
      cp,
      restrict_to_season = FALSE,
      keep_out_of_season_only = FALSE
    ),
    non_rsv_season = clean_and_restrict(
      cp,
      restrict_to_season = FALSE,
      keep_out_of_season_only = TRUE
    )
  )
}

model_data_washout <- map(cp_washout, make_period_datasets)

walk(model_data_washout, function(period_list) {
  walk(period_list, function(md) {
    dup <- md %>%
      filter(event == 1L) %>%
      count(infant_study_id, tstop, name = "n") %>%
      filter(n > 1L)
    stopifnot(nrow(dup) == 0L)
  })
})

fit_one_period <- function(md, period_label, washout_days) {
  fit <- fit_ag(md)

  ve <- tidy_ve(
    fit,
    label = paste0("all_cause_", period_label, "_", washout_days, "d")
  ) %>%
    mutate(
      washout_days = washout_days,
      period = period_label,
      .before = 1
    )

  list(
    fit = fit,
    ve = ve
  )
}

washout_fits <- imap(model_data_washout, function(period_list, washout_name) {
  w <- as.integer(stringr::str_extract(washout_name, "\\d+"))

  imap(period_list, function(md, period_name) {
    fit_one_period(md, period_name, w)
  })
})

washout_ve_results <- imap_dfr(washout_fits, function(period_fits, washout_name) {
  map_dfr(period_fits, "ve")
}) %>%
  mutate(
    period = factor(
      period,
      levels = c("rsv_season", "full_year", "non_rsv_season"),
      labels = c("RSV seasons", "Full year", "Non-RSV seasons")
    ),
    exposure_label = case_when(
      exposure == "exposure_arma_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure == "exposure_armb_nirsevimab_only" ~ "Nirsevimab",
      TRUE ~ exposure
    )
  ) %>%
  arrange(washout_days, period, exposure_label)
print(
  washout_ve_results %>%
    select(
      washout_days, period, exposure_label,
      HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high, p.value
    ),
  n = Inf
)

period_pretty <- c(
  rsv_season = "RSV seasons",
  full_year = "Full year",
  non_rsv_season = "Non-RSV seasons"
)

washout_arm_counts <- imap_dfr(model_data_washout, function(period_list, washout_name) {
  w <- as.integer(stringr::str_extract(washout_name, "\\d+"))

  imap_dfr(period_list, function(md, period_name) {

    check_arm_counts(
      md,
      paste0(w, "-day washout: ", period_pretty[[period_name]])
    )

    arm_counts(md, paste0(w, "-day washout: ", period_pretty[[period_name]])) %>%
      mutate(
        washout_days = w,
        period = period_pretty[[period_name]],
        .before = 1
      )
  })
})
print(
  washout_arm_counts %>%
    select(
      washout_days, period, exposure_label,
      n_infants_contributing, person_days, person_years,
      n_events, events_per_100_py
    ),
  n = Inf
)

stopifnot(
  exists("primary_ve_table"),
  exists("sens_ve_table_fullyear"),
  exists("out_of_season_ve_table"),
  exists("model_data_primary"),
  exists("model_data_fullyear"),
  exists("model_data_out_of_season")
)

original_ve <- bind_rows(
  primary_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "RSV seasons"),
  sens_ve_table_fullyear %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "Full year"),
  out_of_season_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "Non-RSV seasons")
) %>%
  mutate(
    washout_days = 0L,
    definition = "Primary: unique antibiotic-use days",
    exposure_label = case_when(
      exposure == "exposure_arma_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure == "exposure_armb_nirsevimab_only" ~ "Nirsevimab",
      TRUE ~ exposure
    )
  ) %>%
  select(
    washout_days, definition, period, exposure, exposure_label,
    HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high, p.value
  )

washout_ve_for_compare <- washout_ve_results %>%
  mutate(
    period = as.character(period),
    definition = paste0(washout_days, "-day washout")
  ) %>%
  select(
    washout_days, definition, period, exposure, exposure_label,
    HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high, p.value
  )

reviewer_ve_comparison <- bind_rows(
  original_ve,
  washout_ve_for_compare
) %>%
  mutate(
    period = factor(
      period,
      levels = c("RSV seasons", "Full year", "Non-RSV seasons")
    ),
    definition = factor(
      definition,
      levels = c(
        "Primary: unique antibiotic-use days",
        "14-day washout",
        "30-day washout"
      )
    ),
    `HR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    `VE, % (95% CI)` = sprintf(
      "%.1f (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high
    )
  ) %>%
  arrange(period, exposure_label, washout_days)
print(
  reviewer_ve_comparison %>%
    select(
      period, exposure_label, definition,
      `HR (95% CI)`, `VE, % (95% CI)`, p.value
    ),
  n = Inf
)

original_model_sets <- list(
  `RSV seasons` = model_data_primary$all_cause,
  `Full year` = model_data_fullyear$all_cause,
  `Non-RSV seasons` = model_data_out_of_season$all_cause
)

original_arm_counts <- imap_dfr(original_model_sets, function(md, period_name) {
  arm_counts(md, paste0("Primary definition: ", period_name)) %>%
    mutate(
      washout_days = 0L,
      definition = "Primary: unique antibiotic-use days",
      period = period_name,
      .before = 1
    )
})

washout_arm_counts_for_compare <- washout_arm_counts %>%
  mutate(definition = paste0(washout_days, "-day washout"))

reviewer_event_count_comparison <- bind_rows(
  original_arm_counts,
  washout_arm_counts_for_compare
) %>%
  arrange(period, exposure_label, washout_days)
print(
  reviewer_event_count_comparison %>%
    select(
      period, exposure_label, definition,
      n_infants_contributing, person_days, person_years,
      n_events, events_per_100_py
    ),
  n = Inf
)
iwalk(washout_fits, function(period_fits, washout_name) {
  iwalk(period_fits, function(obj, period_name) {
    print(cox.zph(obj$fit))
  })
})

write.csv(
  event_definition_counts,
  "washout_14_30_event_definition_counts.csv",
  row.names = FALSE
)

write.csv(
  washout_arm_counts,
  "washout_14_30_model_specific_event_counts.csv",
  row.names = FALSE
)

write.csv(
  washout_ve_results,
  "washout_14_30_adjusted_VE_results.csv",
  row.names = FALSE
)

write.csv(
  reviewer_ve_comparison,
  "washout_primary_14_30_VE_comparison.csv",
  row.names = FALSE
)

write.csv(
  reviewer_event_count_comparison,
  "washout_primary_14_30_event_count_comparison.csv",
  row.names = FALSE
)

saveRDS(
  washout_fits,
  "washout_14_30_fitted_models.rds"
)
washout_counts <- read.csv(
  "washout_14_30_model_specific_event_counts.csv"
)

names(washout_counts)

print(washout_counts, row.names = FALSE)

library(dplyr)
library(tidyr)

original_events <- events_primary$all_cause %>%
  distinct(infant_study_id, abx_days) %>%
  arrange(infant_study_id, abx_days)

events_14 <- apply_washout(original_events, 14)
events_30 <- apply_washout(original_events, 30)

collapsed_14 <- original_events %>%
  anti_join(
    events_14,
    by = c("infant_study_id", "abx_days")
  ) %>%
  mutate(washout_days = 14)

collapsed_30 <- original_events %>%
  anti_join(
    events_30,
    by = c("infant_study_id", "abx_days")
  ) %>%
  mutate(washout_days = 30)

collapsed_events <- bind_rows(
  collapsed_14,
  collapsed_30
)
collapsed_events %>%
  count(washout_days) %>%
  print()

infant_dates_full <- cohort_step2 %>%
  select(infant_study_id, index_date) %>%
  distinct()

collapsed_events_dated <- collapsed_events %>%
  left_join(
    infant_dates_full,
    by = "infant_study_id"
  ) %>%
  mutate(
    event_date = index_date + abx_days
  )
summary(collapsed_events_dated$event_date)

collapsed_events_dated <- collapsed_events_dated %>%
  mutate(
    period = case_when(
      event_date >= as.Date("2023-10-01") &
        event_date <= as.Date("2024-02-29") ~ "RSV season",

      event_date >= as.Date("2024-11-01") &
        event_date <= as.Date("2025-03-31") ~ "RSV season",

      TRUE ~ "Non-RSV season"
    )
  )

collapsed_by_period <- collapsed_events_dated %>%
  count(washout_days, period) %>%
  group_by(washout_days) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(collapsed_by_period, n = Inf)

season_boundaries <- as.Date(c(
  "2023-10-01",
  "2024-03-01",
  "2024-11-01",
  "2025-04-01"
))

collapsed_boundary <- collapsed_events_dated %>%
  rowwise() %>%
  mutate(
    days_from_nearest_boundary =
      min(abs(as.numeric(event_date - season_boundaries)))
  ) %>%
  ungroup() %>%
  mutate(
    boundary_window = case_when(
      days_from_nearest_boundary <= 14 ~ "0-14 days",
      days_from_nearest_boundary <= 30 ~ "15-30 days",
      TRUE ~ ">30 days"
    )
  )

boundary_summary <- collapsed_boundary %>%
  count(washout_days, boundary_window) %>%
  group_by(washout_days) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(boundary_summary, n = Inf)

event_history <- original_events %>%
  group_by(infant_study_id) %>%
  arrange(abx_days, .by_group = TRUE) %>%
  mutate(
    previous_abx_day = lag(abx_days),
    gap_days = abx_days - previous_abx_day
  ) %>%
  ungroup() %>%
  left_join(
    infant_dates_full,
    by = "infant_study_id"
  ) %>%
  mutate(
    event_date = index_date + abx_days,
    previous_event_date = index_date + previous_abx_day,

    event_period = case_when(
      event_date >= as.Date("2023-10-01") &
        event_date <= as.Date("2024-02-29") ~ "RSV season",

      event_date >= as.Date("2024-11-01") &
        event_date <= as.Date("2025-03-31") ~ "RSV season",

      TRUE ~ "Non-RSV season"
    ),

    previous_period = case_when(
      is.na(previous_event_date) ~ NA_character_,

      previous_event_date >= as.Date("2023-10-01") &
        previous_event_date <= as.Date("2024-02-29") ~ "RSV season",

      previous_event_date >= as.Date("2024-11-01") &
        previous_event_date <= as.Date("2025-03-31") ~ "RSV season",

      TRUE ~ "Non-RSV season"
    )
  )

cross_boundary_events <- bind_rows(

  event_history %>%
    filter(!is.na(gap_days), gap_days <= 14) %>%
    mutate(washout_days = 14),

  event_history %>%
    filter(!is.na(gap_days), gap_days <= 30) %>%
    mutate(washout_days = 30)

) %>%
  mutate(
    crosses_period = event_period != previous_period
  )

cross_boundary_summary <- cross_boundary_events %>%
  count(washout_days, crosses_period) %>%
  group_by(washout_days) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(cross_boundary_summary, n = Inf)

library(dplyr)
library(tidyr)

recurrence_events <- events_primary$all_cause %>%
  distinct(infant_study_id, abx_days) %>%
  arrange(infant_study_id, abx_days) %>%
  group_by(infant_study_id) %>%
  mutate(
    previous_abx_day = lag(abx_days),
    gap_days = abx_days - previous_abx_day
  ) %>%
  ungroup() %>%
  left_join(
    infant_dates_full,
    by = "infant_study_id"
  ) %>%
  mutate(
    event_date = index_date + abx_days,

    period = case_when(
      event_date >= as.Date("2023-10-01") &
        event_date <= as.Date("2024-02-29") ~ "RSV seasons",

      event_date >= as.Date("2024-11-01") &
        event_date <= as.Date("2025-03-31") ~ "RSV seasons",

      TRUE ~ "Non-RSV seasons"
    ),

    gap_category = case_when(
      is.na(gap_days) ~ "First event",
      gap_days <= 14 ~ "1-14 days",
      gap_days <= 30 ~ "15-30 days",
      TRUE ~ ">30 days"
    )
  )

recurrence_events %>%
  count(gap_category) %>%
  print(n = Inf)

names(cp_primary)

lapply(cp_primary, function(x) {
  c(
    rows = nrow(x),
    infants = n_distinct(x$infant_study_id),
    events = sum(x$event)
  )
})

recurrence_events_exposed <- recurrence_events %>%
  inner_join(
    event_exposure_lookup,
    by = c("infant_study_id", "abx_days")
  ) %>%
  mutate(
    exposure_label = case_when(
      exposure_arm == "a_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure_arm == "b_nirsevimab_only" ~ "Nirsevimab",
      exposure_arm == "d_neither" ~ "Neither",
      TRUE ~ exposure_arm
    )
  )

recurrence_events_exposed %>%
  count(exposure_label) %>%
  print(n = Inf)

recurrence_pattern <- recurrence_events_exposed %>%
  count(
    period,
    exposure_label,
    gap_category,
    name = "n_events"
  ) %>%
  group_by(
    period,
    exposure_label
  ) %>%
  mutate(
    total_events = sum(n_events),
    percent = 100 * n_events / total_events
  ) %>%
  ungroup() %>%
  arrange(
    period,
    exposure_label,
    gap_category
  )

print(recurrence_pattern, n = Inf)

recurrence_summary <- recurrence_events_exposed %>%
  group_by(
    period,
    exposure_label
  ) %>%
  summarise(
    total_antibiotic_days = n(),

    first_events =
      sum(is.na(gap_days)),

    repeat_within_14 =
      sum(!is.na(gap_days) & gap_days <= 14),

    repeat_within_30 =
      sum(!is.na(gap_days) & gap_days <= 30),

    pct_repeat_within_14 =
      100 * repeat_within_14 / total_antibiotic_days,

    pct_repeat_within_30 =
      100 * repeat_within_30 / total_antibiotic_days,

    .groups = "drop"
  ) %>%
  arrange(
    period,
    exposure_label
  )

print(recurrence_summary, n = Inf)

write.csv(
  recurrence_summary,
  "washout_recurrence_pattern_summary.csv",
  row.names = FALSE
)

write.csv(
  recurrence_pattern,
  "washout_recurrence_gap_distribution.csv",
  row.names = FALSE
)
