

library(tidyverse)
library(survival)
library(broom)

stopifnot(
  exists("analytic_events"),
  exists("cohort_primary"),
  exists("FOLLOWUP_MAX"),
  exists("blocks_primary"),
  exists("build_cp"),
  exists("clean_and_restrict"),
  exists("fit_ag"),
  exists("tidy_ve"),
  exists("arm_counts"),
  exists("COVARIATE_TERMS")
)

EXCLUDED_ABX_CLASSES <- c(
  "Aminoglycoside Antibiotics",
  "Vancomycin Antibiotics And Derivatives"
)
available_classes <- analytic_events %>%
  filter(!is.na(abx_class), abx_class != "") %>%
  distinct(abx_class) %>%
  pull(abx_class) %>%
  as.character()

missing_classes <- setdiff(EXCLUDED_ABX_CLASSES, available_classes)

stopifnot(
  "One or more excluded antibiotic classes were not found in analytic_events$abx_class" =
    length(missing_classes) == 0
)
print(EXCLUDED_ABX_CLASSES)

eligible_abx_records <- analytic_events %>%
  filter(
    analyzable_event_primary == TRUE,
    !is.na(abx_days),
    abx_days > 0,
    abx_days <= FOLLOWUP_MAX
  ) %>%
  semi_join(
    cohort_primary %>% distinct(infant_study_id),
    by = "infant_study_id"
  )

all_cause_event_days_reconstructed <- eligible_abx_records %>%
  distinct(infant_study_id, abx_days)

all_cause_no_amino_vanc_events <- eligible_abx_records %>%
  filter(
    is.na(abx_class) | !(abx_class %in% EXCLUDED_ABX_CLASSES)
  ) %>%
  distinct(infant_study_id, abx_days)

if (exists("events_primary") &&
    is.list(events_primary) &&
    "all_cause" %in% names(events_primary)) {

  primary_event_days <- events_primary$all_cause %>%
    filter(
      !is.na(abx_days),
      abx_days > 0,
      abx_days <= FOLLOWUP_MAX
    ) %>%
    semi_join(
      cohort_primary %>% distinct(infant_study_id),
      by = "infant_study_id"
    ) %>%
    distinct(infant_study_id, abx_days)

  validation_primary <- list(
    n_primary_object = nrow(primary_event_days),
    n_reconstructed  = nrow(all_cause_event_days_reconstructed),
    only_in_primary  = nrow(anti_join(
      primary_event_days,
      all_cause_event_days_reconstructed,
      by = c("infant_study_id", "abx_days")
    )),
    only_in_reconstructed = nrow(anti_join(
      all_cause_event_days_reconstructed,
      primary_event_days,
      by = c("infant_study_id", "abx_days")
    ))
  )
  print(validation_primary)

  stopifnot(
    "Reconstructed all-cause event days do not match events_primary$all_cause" =
      validation_primary$only_in_primary == 0 &&
      validation_primary$only_in_reconstructed == 0
  )
} else {
  warning(
    "events_primary$all_cause was not available; reconstruction could not be ",
    "cross-checked against the primary event object."
  )
}

removed_event_days <- anti_join(
  all_cause_event_days_reconstructed,
  all_cause_no_amino_vanc_events,
  by = c("infant_study_id", "abx_days")
)

days_with_excluded_class <- eligible_abx_records %>%
  filter(abx_class %in% EXCLUDED_ABX_CLASSES) %>%
  distinct(infant_study_id, abx_days)

mixed_event_days <- inner_join(
  days_with_excluded_class,
  all_cause_no_amino_vanc_events,
  by = c("infant_study_id", "abx_days")
)

event_exclusion_summary <- tibble(
  metric = c(
    "Primary all-cause infant-days",
    "Infant-days containing aminoglycoside and/or vancomycin",
    "Mixed days also containing another retained antibiotic",
    "All-cause infant-days removed after exclusion",
    "All-cause infant-days retained"
  ),
  n = c(
    nrow(all_cause_event_days_reconstructed),
    nrow(days_with_excluded_class),
    nrow(mixed_event_days),
    nrow(removed_event_days),
    nrow(all_cause_no_amino_vanc_events)
  )
)
print(event_exclusion_summary, n = Inf)
cp_no_amino_vanc <- build_cp(
  all_cause_no_amino_vanc_events,
  blocks = blocks_primary
)

md_no_amino_vanc <- list(
  RSV_season = clean_and_restrict(
    cp_no_amino_vanc,
    restrict_to_season = TRUE
  ),
  Full_year = clean_and_restrict(
    cp_no_amino_vanc,
    restrict_to_season = FALSE
  ),
  Non_RSV_season = clean_and_restrict(
    cp_no_amino_vanc,
    keep_out_of_season_only = TRUE
  )
)

walk(md_no_amino_vanc, function(md) {
  stopifnot(
    all(md$tstop > md$tstart),
    all(md$tstart >= 0),
    all(md$tstop <= FOLLOWUP_MAX)
  )
})

fit_sensitivity_model <- function(md, period_label) {

  fit <- fit_ag(md)

  ve <- tidy_ve(
    fit,
    paste0("All-cause excluding aminoglycoside/vancomycin - ", period_label)
  )

  counts <- arm_counts(
    md,
    paste0("All-cause excluding aminoglycoside/vancomycin - ", period_label)
  )

  list(
    fit = fit,
    ve = ve,
    counts = counts
  )
}

sensitivity_fits <- imap(
  md_no_amino_vanc,
  ~ fit_sensitivity_model(.x, .y)
)

model_results_no_amino_vanc <- imap_dfr(
  sensitivity_fits,
  function(x, period) {
    x$ve %>%
      mutate(analysis_period = period, .before = 1)
  }
)
print(
  model_results_no_amino_vanc %>%
    select(
      any_of(c(
        "analysis_period", "outcome", "exposure",
        "HR", "CI_low", "CI_high",
        "VE_pct", "VE_CI_low", "VE_CI_high",
        "p.value"
      ))
    ),
  n = Inf
)

arm_counts_no_amino_vanc <- imap_dfr(
  sensitivity_fits,
  function(x, period) {
    x$counts %>%
      mutate(analysis_period = period, .before = 1)
  }
)
print(arm_counts_no_amino_vanc, n = Inf)

overall_counts_no_amino_vanc <- imap_dfr(
  md_no_amino_vanc,
  function(md, period) {
    tibble(
      analysis_period = period,
      n_infants = n_distinct(md$infant_study_id),
      n_mothers = n_distinct(md$mom_study_id),
      n_intervals = nrow(md),
      person_days = sum(md$tstop - md$tstart),
      n_events = sum(md$event)
    )
  }
)
print(overall_counts_no_amino_vanc, n = Inf)

if (exists("cp_primary") &&
    is.list(cp_primary) &&
    "all_cause" %in% names(cp_primary)) {

  md_original_all_cause <- list(
    RSV_season = clean_and_restrict(
      cp_primary$all_cause,
      restrict_to_season = TRUE
    ),
    Full_year = clean_and_restrict(
      cp_primary$all_cause,
      restrict_to_season = FALSE
    ),
    Non_RSV_season = clean_and_restrict(
      cp_primary$all_cause,
      keep_out_of_season_only = TRUE
    )
  )

  original_results <- imap_dfr(
    md_original_all_cause,
    function(md, period) {
      fit <- fit_ag(md)
      tidy_ve(fit, paste0("Original all-cause - ", period)) %>%
        mutate(
          analysis_period = period,
          analysis = "Original all-cause",
          .before = 1
        )
    }
  )

  excluded_results <- model_results_no_amino_vanc %>%
    mutate(
      analysis = "Excluding aminoglycoside/vancomycin",
      .after = analysis_period
    )

  comparison_results <- bind_rows(
    original_results,
    excluded_results
  ) %>%
    mutate(
      exposure_label = case_when(
        str_detect(exposure, "maternal_vax_only") ~ "Maternal RSV vaccine",
        str_detect(exposure, "nirsevimab_only")   ~ "Nirsevimab",
        TRUE                                      ~ exposure
      )
    ) %>%
    select(
      analysis_period,
      analysis,
      exposure_label,
      any_of(c(
        "HR", "CI_low", "CI_high",
        "VE_pct", "VE_CI_low", "VE_CI_high",
        "p.value"
      ))
    ) %>%
    arrange(
      factor(analysis_period,
             levels = c("Full_year", "RSV_season", "Non_RSV_season")),
      exposure_label,
      analysis
    )
  print(comparison_results, n = Inf)

} else {
  comparison_results <- NULL
  warning(
    "cp_primary$all_cause not available; original-vs-sensitivity comparison ",
    "was not generated."
  )
}

write_csv(
  event_exclusion_summary,
  "exclude_amino_vanc_event_summary.csv"
)

write_csv(
  model_results_no_amino_vanc,
  "exclude_amino_vanc_model_results.csv"
)

write_csv(
  arm_counts_no_amino_vanc,
  "exclude_amino_vanc_arm_counts.csv"
)

write_csv(
  overall_counts_no_amino_vanc,
  "exclude_amino_vanc_overall_counts.csv"
)

if (!is.null(comparison_results)) {
  write_csv(
    comparison_results,
    "exclude_amino_vanc_comparison.csv"
  )
}
if (!is.null(comparison_results)) {
}
