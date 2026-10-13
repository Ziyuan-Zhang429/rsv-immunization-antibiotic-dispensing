

rm(list = ls())

PRIMARY_SCRIPT <- "R/primary/cox_primary.R"

if (!file.exists(PRIMARY_SCRIPT)) {
  stop(
    "Cannot find ", PRIMARY_SCRIPT,
    ". Put this script in the same working directory as the primary analysis ",
    "script and SAS input files."
  )
}

source(PRIMARY_SCRIPT)

required_objects <- c(
  "cohort_primary", "eligible_cohort", "events_primary",
  "make_analysis_blocks", "build_cp", "clean_and_restrict",
  "fit_ag", "tidy_ve"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "Primary script did not create required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

all_cause_events <- events_primary$all_cause

fit_three_periods <- function(cp, sensitivity_label) {

  period_data <- list(
    "Full year" = clean_and_restrict(
      cp,
      restrict_to_season = FALSE
    ),

    "RSV season" = clean_and_restrict(
      cp,
      restrict_to_season = TRUE
    ),

    "Non-RSV season" = clean_and_restrict(
      cp,
      keep_out_of_season_only = TRUE
    )
  )

  results <- purrr::imap_dfr(
    period_data,
    function(md, period_name) {
      fit <- fit_ag(md)

      tidy_ve(fit, "all_cause") %>%
        mutate(
          sensitivity = sensitivity_label,
          period = period_name,
          .before = 1
        )
    }
  )

  list(
    data = period_data,
    results = results
  )
}

pretty_exposure <- function(x) {
  case_when(
    grepl("a_maternal_vax_only", x) ~ "Maternal RSVpreF",
    grepl("b_nirsevimab_only", x)    ~ "Nirsevimab",
    TRUE                             ~ x
  )
}

PRODUCT_AVAIL_BIRTH_CUTOFF <- as.Date("2023-11-01")
cohort_availability <- cohort_primary %>%
  filter(index_date >= PRODUCT_AVAIL_BIRTH_CUTOFF)
blocks_availability <- make_analysis_blocks(
  cohort_availability,
  label = "R2_product_availability"
)

cp_availability <- build_cp(
  all_cause_events,
  blocks = blocks_availability,
  maternal_exposure_flag = "rsv_valid"
)

availability_fit <- fit_three_periods(
  cp_availability,
  "Births >= Nov 1, 2023"
)

results_availability <- availability_fit$results
print(results_availability, n = Inf)
cohort_anytime_all <- eligible_cohort %>%
  mutate(
    serious_newborn_comorbidity = as.integer(
      bronchopulm_dx == 1 | immunodef_dx == 1 | downsyn_dx == 1
    ),

    rsv_valid_anytime =
      rsv_history &
      !sibling_match_confirmed &
      RSV_days > 0,

    ever_both_anytime =
      rsv_valid_anytime & mab_history
  )

maternal_reclassification <- cohort_anytime_all %>%
  count(rsv_valid, rsv_valid_anytime, name = "n") %>%
  mutate(pct = 100 * n / sum(n))
print(maternal_reclassification, n = Inf)

newly_exposed_anytime <- cohort_anytime_all %>%
  filter(!rsv_valid, rsv_valid_anytime)
newly_exposed_timing_summary <- newly_exposed_anytime %>%
  summarise(
    n = n(),

    n_vax_lt_14d_before_birth =
      sum(
        RSV_days > 0 &
          RSV_days < MIN_DAYS_VAX_TO_BIRTH,
        na.rm = TRUE
      ),

    n_vax_before_32wk =
      sum(
        gestage_at_vax < GESTAGE_AT_VAX_MIN,
        na.rm = TRUE
      ),

    n_vax_after_36wk =
      sum(
        gestage_at_vax > GESTAGE_AT_VAX_MAX,
        na.rm = TRUE
      ),

    median_days_vax_to_birth =
      median(RSV_days, na.rm = TRUE),

    min_days_vax_to_birth =
      min(RSV_days, na.rm = TRUE),

    max_days_vax_to_birth =
      max(RSV_days, na.rm = TRUE)
  )
print(newly_exposed_timing_summary)

cohort_anytime <- cohort_anytime_all %>%
  filter(!ever_both_anytime)
blocks_anytime <- make_analysis_blocks(
  cohort_anytime,
  label = "R2_maternal_anytime"
)

cp_anytime <- build_cp(
  all_cause_events,
  blocks = blocks_anytime,
  maternal_exposure_flag = "rsv_valid_anytime"
)

timing_fit <- fit_three_periods(
  cp_anytime,
  "Any maternal RSVpreF before delivery"
)

results_timing <- timing_fit$results
print(results_timing, n = Inf)

reviewer2_sensitivity_results <- bind_rows(
  results_availability,
  results_timing
) %>%
  mutate(
    exposure = pretty_exposure(exposure),
    period = factor(
      period,
      levels = c("Full year", "RSV season", "Non-RSV season")
    )
  ) %>%
  arrange(sensitivity, period, exposure) %>%
  select(
    sensitivity,
    period,
    exposure,
    VE_pct,
    VE_CI_low,
    VE_CI_high,
    HR,
    CI_low,
    CI_high,
    p.value
  )
print(reviewer2_sensitivity_results, n = Inf)

primary_allcause_comparison <- bind_rows(
  sens_ve_table_fullyear %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "Full year"),

  primary_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "RSV season"),

  out_of_season_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period = "Non-RSV season")
) %>%
  mutate(
    sensitivity = "Primary analysis",
    exposure = pretty_exposure(exposure)
  ) %>%
  select(
    sensitivity,
    period,
    exposure,
    VE_pct,
    VE_CI_low,
    VE_CI_high,
    HR,
    CI_low,
    CI_high,
    p.value
  )

reviewer2_with_primary <- bind_rows(
  primary_allcause_comparison,
  reviewer2_sensitivity_results
) %>%
  mutate(
    period = factor(
      period,
      levels = c("Full year", "RSV season", "Non-RSV season")
    )
  ) %>%
  arrange(
    factor(
      sensitivity,
      levels = c(
        "Primary analysis",
        "Births >= Nov 1, 2023",
        "Any maternal RSVpreF before delivery"
      )
    ),
    period,
    exposure
  )
print(reviewer2_with_primary, n = Inf)

write.csv(
  reviewer2_sensitivity_results,
  "Reviewer2_immortal_time_sensitivity_results.csv",
  row.names = FALSE
)

write.csv(
  reviewer2_with_primary,
  "Reviewer2_primary_and_sensitivity_comparison.csv",
  row.names = FALSE
)

write.csv(
  maternal_reclassification,
  "Reviewer2_maternal_RSVpreF_reclassification.csv",
  row.names = FALSE
)

write.csv(
  newly_exposed_timing_summary,
  "Reviewer2_maternal_RSVpreF_newly_exposed_timing_summary.csv",
  row.names = FALSE
)
