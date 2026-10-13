

stopifnot(
  exists("cp_primary"),
  exists("clean_and_restrict"),
  exists("fit_ag"),
  exists("tidy_ve"),
  exists("arm_counts"),
  exists("check_arm_counts"),
  exists("FOLLOWUP_MAX")
)

library(tidyverse)

EARLY_EXCLUSION_WINDOWS <- tribble(
  ~window_code,       ~cutoff, ~analysis_window,
  "primary",                0, "Days 1-180",
  "exclude_days_1_7",       7, "Days 8-180",
  "exclude_days_1_14",     14, "Days 15-180",
  "exclude_days_1_28",     28, "Days 29-180"
)

print(EARLY_EXCLUSION_WINDOWS)

apply_delayed_entry <- function(md, cutoff) {
  stopifnot(
    is.numeric(cutoff),
    length(cutoff) == 1,
    cutoff >= 0,
    cutoff < FOLLOWUP_MAX
  )

  if (cutoff == 0) return(md)

  out <- md %>%

    filter(tstop > cutoff) %>%

    mutate(
      tstart = pmax(tstart, cutoff),

      calendar_date_start = index_date + tstart,
      calendar_month = floor_date(calendar_date_start, "month")
    ) %>%

    filter(tstop > tstart) %>%
    mutate(
      calendar_month = droplevels(factor(calendar_month))
    )

  stopifnot(
    all(out$tstart >= cutoff),
    all(out$tstop <= FOLLOWUP_MAX),
    all(out$tstop > out$tstart)
  )

  out
}

early_md_rsv <- map(
  EARLY_EXCLUSION_WINDOWS$cutoff,
  ~ cp_primary$all_cause %>%
    clean_and_restrict(
      restrict_to_season = TRUE,
      keep_out_of_season_only = FALSE
    ) %>%
    apply_delayed_entry(.x)
)
names(early_md_rsv) <- EARLY_EXCLUSION_WINDOWS$window_code

early_md_full <- map(
  EARLY_EXCLUSION_WINDOWS$cutoff,
  ~ cp_primary$all_cause %>%
    clean_and_restrict(
      restrict_to_season = FALSE,
      keep_out_of_season_only = FALSE
    ) %>%
    apply_delayed_entry(.x)
)
names(early_md_full) <- EARLY_EXCLUSION_WINDOWS$window_code

early_md_oos <- map(
  EARLY_EXCLUSION_WINDOWS$cutoff,
  ~ cp_primary$all_cause %>%
    clean_and_restrict(
      restrict_to_season = FALSE,
      keep_out_of_season_only = TRUE
    ) %>%
    apply_delayed_entry(.x)
)
names(early_md_oos) <- EARLY_EXCLUSION_WINDOWS$window_code

validate_delayed_entry <- function(md_list, period_label) {
  imap(md_list, function(md, nm) {
    cutoff <- EARLY_EXCLUSION_WINDOWS %>%
      filter(window_code == nm) %>%
      pull(cutoff)

    stopifnot(
      all(md$tstart >= cutoff),
      all(md$tstop <= FOLLOWUP_MAX),
      all(md$tstop > md$tstart)
    )

    duplicated_events <- md %>%
      filter(event == 1L) %>%
      count(infant_study_id, tstop, name = "n") %>%
      filter(n > 1)

    if (nrow(duplicated_events) > 0) {
      stop(period_label, " / ", nm,
           ": multiple antibiotic events found for the same infant-day.")
    }

    invisible(TRUE)
  })
}

validate_delayed_entry(early_md_rsv,  "RSV seasons")
validate_delayed_entry(early_md_full, "Full year")
validate_delayed_entry(early_md_oos,  "Non-RSV seasons")

summarize_window <- function(md_list, period_label) {
  imap_dfr(md_list, function(md, nm) {
    spec <- EARLY_EXCLUSION_WINDOWS %>% filter(window_code == nm)

    tibble(
      period = period_label,
      window_code = nm,
      analysis_window = spec$analysis_window,
      cutoff = spec$cutoff,
      n_infants = n_distinct(md$infant_study_id),
      person_days = sum(md$tstop - md$tstart),
      n_events = sum(md$event)
    )
  })
}

early_window_overall_counts <- bind_rows(
  summarize_window(early_md_rsv,  "RSV seasons"),
  summarize_window(early_md_full, "Full year"),
  summarize_window(early_md_oos,  "Non-RSV seasons")
)
print(early_window_overall_counts, n = Inf)

get_arm_counts_windows <- function(md_list, period_label) {
  imap_dfr(md_list, function(md, nm) {
    spec <- EARLY_EXCLUSION_WINDOWS %>% filter(window_code == nm)

    invisible(
      check_arm_counts(
        md,
        paste(period_label, spec$analysis_window, sep = " | ")
      )
    )

    arm_counts(
      md,
      paste(period_label, spec$analysis_window, sep = " | ")
    ) %>%
      mutate(
        period = period_label,
        window_code = nm,
        analysis_window = spec$analysis_window,
        cutoff = spec$cutoff,
        .before = 1
      )
  })
}

early_counts <- bind_rows(
  get_arm_counts_windows(early_md_rsv,  "RSV seasons"),
  get_arm_counts_windows(early_md_full, "Full year"),
  get_arm_counts_windows(early_md_oos,  "Non-RSV seasons")
)
print(
  early_counts %>%
    select(period, analysis_window, exposure_label,
           n_infants_contributing, person_days, person_years, n_events),
  n = Inf
)

fit_window_models <- function(md_list, period_label) {
  imap_dfr(md_list, function(md, nm) {
    spec <- EARLY_EXCLUSION_WINDOWS %>% filter(window_code == nm)

    fit <- fit_ag(md)

    tidy_ve(fit, label = "all_cause") %>%
      mutate(
        period = period_label,
        window_code = nm,
        analysis_window = spec$analysis_window,
        cutoff = spec$cutoff,
        .before = 1
      )
  })
}

early_ve_results <- bind_rows(
  fit_window_models(early_md_rsv,  "RSV seasons"),
  fit_window_models(early_md_full, "Full year"),
  fit_window_models(early_md_oos,  "Non-RSV seasons")
)
print(
  early_ve_results %>%
    select(period, analysis_window, exposure, HR, CI_low, CI_high,
           VE_pct, VE_CI_low, VE_CI_high, p.value),
  n = Inf
)

early_ve_exposed <- early_ve_results %>%
  filter(
    exposure %in% c(
      "exposure_arma_maternal_vax_only",
      "exposure_armb_nirsevimab_only"
    )
  ) %>%
  mutate(
    exposure_label = case_when(
      exposure == "exposure_arma_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure == "exposure_armb_nirsevimab_only"   ~ "Nirsevimab"
    )
  )

early_sensitivity_table <- early_ve_exposed %>%
  left_join(
    early_counts %>%
      select(period, window_code, exposure, n_infants_contributing,
             person_days, person_years, n_events),
    by = c("period", "window_code", "exposure")
  ) %>%
  mutate(
    `Analysis period` = period,
    `Follow-up window` = analysis_window,
    Exposure = exposure_label,
    `HR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high),
    `VE, % (95% CI)` = sprintf("%.1f (%.1f, %.1f)",
                               VE_pct, VE_CI_low, VE_CI_high)
  ) %>%
  select(
    `Analysis period`,
    `Follow-up window`,
    Exposure,
    `Infants contributing person-time` = n_infants_contributing,
    `Person-days` = person_days,
    `Antibiotic events` = n_events,
    `HR (95% CI)`,
    `VE, % (95% CI)`
  ) %>%
  arrange(
    factor(`Analysis period`,
           levels = c("Full year", "RSV seasons", "Non-RSV seasons")),
    factor(`Follow-up window`,
           levels = EARLY_EXCLUSION_WINDOWS$analysis_window),
    factor(Exposure,
           levels = c("Maternal RSV vaccine", "Nirsevimab"))
  )
print(early_sensitivity_table, n = Inf)

early_ve_compact <- early_ve_exposed %>%
  mutate(
    Exposure = exposure_label,
    `VE, % (95% CI)` = sprintf("%.1f (%.1f, %.1f)",
                               VE_pct, VE_CI_low, VE_CI_high)
  ) %>%
  select(period, analysis_window, Exposure, `VE, % (95% CI)`) %>%
  pivot_wider(
    names_from = Exposure,
    values_from = `VE, % (95% CI)`
  ) %>%
  arrange(
    factor(period,
           levels = c("Full year", "RSV seasons", "Non-RSV seasons")),
    factor(analysis_window,
           levels = EARLY_EXCLUSION_WINDOWS$analysis_window)
  )
print(early_ve_compact, n = Inf)

event_loss_table <- early_window_overall_counts %>%
  group_by(period) %>%
  mutate(
    primary_events = n_events[window_code == "primary"],
    events_removed = primary_events - n_events,
    pct_events_removed = 100 * events_removed / primary_events
  ) %>%
  ungroup() %>%
  select(period, analysis_window, n_infants, person_days, n_events,
         events_removed, pct_events_removed)
print(event_loss_table, n = Inf)

early_models_for_ph <- list(
  rsv  = early_md_rsv,
  full = early_md_full,
  oos  = early_md_oos
)

early_ph_results <- imap(
  early_models_for_ph,
  function(md_list, period_code) {
    imap(
      md_list,
      function(md, window_code) {
        fit <- fit_ag(md)
        list(
          fit = fit,
          zph = cox.zph(fit)
        )
      }
    )
  }
)
walk(names(early_ph_results), function(period_code) {
  walk(names(early_ph_results[[period_code]]), function(window_code) {
    print(early_ph_results[[period_code]][[window_code]]$zph)
  })
})

write_csv(
  early_sensitivity_table,
  "early_life_exclusion_sensitivity_table.csv"
)

write_csv(
  early_ve_results,
  "early_life_exclusion_model_results.csv"
)

write_csv(
  early_window_overall_counts,
  "early_life_exclusion_overall_counts.csv"
)

write_csv(
  early_counts,
  "early_life_exclusion_arm_counts.csv"
)

write_csv(
  event_loss_table,
  "early_life_exclusion_event_loss.csv"
)

check_monotonic <- early_window_overall_counts %>%
  arrange(period, cutoff) %>%
  group_by(period) %>%
  summarise(
    person_time_monotonic = all(diff(person_days) <= 0),
    events_monotonic      = all(diff(n_events) <= 0),
    .groups = "drop"
  )
print(check_monotonic)

stopifnot(
  all(check_monotonic$person_time_monotonic),
  all(check_monotonic$events_monotonic)
)
