

library(tidyverse)
library(survival)
library(broom)

stopifnot(
  exists("cp_primary"),
  exists("model_data_primary"),
  exists("model_data_fullyear"),
  exists("model_data_out_of_season"),
  exists("COVARIATE_TERMS"),
  exists("fit_ag"),
  exists("tidy_ve"),
  exists("analytic_events"),
  exists("build_cp"),
  exists("clean_and_restrict")
)

COV_FULL <- COVARIATE_TERMS

COV_NO_MATERNAL_VAX <- paste(
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED +",
  "tot_covid + gbs_prior + utd_status + lb_prior"
)

COV_NO_INFANT_VAX <- paste(
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED + flu_prior + tdap_prior +",
  "tot_covid + gbs_prior + lb_prior"
)

COV_NO_VAX <- paste(
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED +",
  "tot_covid + gbs_prior + lb_prior"
)

make_ag_formula <- function(covariates) {
  as.formula(paste(
    "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_month) +",
    covariates
  ))
}

FORMULAS_R3 <- list(
  `Primary: fully adjusted` = make_ag_formula(COV_FULL),
  `No maternal flu/Tdap` = make_ag_formula(COV_NO_MATERNAL_VAX),
  `No infant routine vaccination` = make_ag_formula(COV_NO_INFANT_VAX),
  `No maternal flu/Tdap or infant vaccination` = make_ag_formula(COV_NO_VAX)
)

restrict_age_window <- function(md, age_start = 0, age_end = 180) {
  md %>%
    filter(tstart >= age_start, tstop <= age_end)
}

R3_WINDOWS <- tribble(
  ~window,                 ~dataset,         ~age_start, ~age_end,
  "Full year, days 1-180", "full_year",      0,          180,
  "RSV season, days 1-180","rsv_season",     0,          180,
  "Non-RSV season, days 1-180", "non_rsv",   0,          180,
  "Full year, days 1-90",  "full_year",      0,           90,
  "RSV season, days 1-90", "rsv_season",     0,           90,
  "Full year, days 91-180", "full_year",    90,          180,
  "RSV season, days 91-180","rsv_season",   90,          180
)

BASE_DATASETS <- list(
  full_year  = model_data_fullyear$all_cause,
  rsv_season = model_data_primary$all_cause,
  non_rsv    = model_data_out_of_season$all_cause
)

fit_r3_ve <- function(md, formula, model_label, window_label) {
  fit <- fit_ag(md, formula = formula)
  tidy_ve(fit, "all_cause") %>%
    mutate(
      window = window_label,
      adjustment = model_label,
      n_infants = n_distinct(md$infant_study_id),
      n_events = sum(md$event),
      .before = 1
    )
}

r3_ve_results <- pmap_dfr(R3_WINDOWS, function(window, dataset, age_start, age_end) {
  md <- restrict_age_window(BASE_DATASETS[[dataset]], age_start, age_end)

  imap_dfr(FORMULAS_R3, function(formula, model_label) {
    fit_r3_ve(md, formula, model_label, window)
  })
}) %>%
  mutate(
    exposure_label = case_when(
      str_detect(exposure, "a_maternal_vax_only") ~ "Maternal RSVpreF",
      str_detect(exposure, "b_nirsevimab_only")   ~ "Nirsevimab",
      TRUE ~ exposure
    ),
    VE_95CI = sprintf("%.1f%% (%.1f to %.1f)", VE_pct, VE_CI_low, VE_CI_high)
  ) %>%
  select(window, exposure_label, adjustment, n_infants, n_events,
         HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high, VE_95CI, p.value)
print(r3_ve_results, n = Inf)

r3_ve_table_wide <- r3_ve_results %>%
  select(window, exposure_label, adjustment, VE_95CI) %>%
  pivot_wider(names_from = adjustment, values_from = VE_95CI) %>%
  arrange(
    factor(window, levels = R3_WINDOWS$window),
    factor(exposure_label, levels = c("Maternal RSVpreF", "Nirsevimab"))
  )
print(r3_ve_table_wide, n = Inf, width = Inf)

r3_ve_change <- r3_ve_results %>%
  select(window, exposure_label, adjustment, VE_pct) %>%
  group_by(window, exposure_label) %>%
  mutate(
    VE_primary = VE_pct[adjustment == "Primary: fully adjusted"],
    delta_VE_pp = VE_pct - VE_primary
  ) %>%
  ungroup() %>%
  arrange(window, exposure_label, adjustment)
print(r3_ve_change, n = Inf)

TARGET_CLASSES_R3 <- c(
  "Aminoglycoside Antibiotics",
  "Vancomycin Antibiotics And Derivatives"
)

class_levels_r3 <- analytic_events %>%
  filter(
    analyzable_event_primary == TRUE,
    !is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX,
    !is.na(abx_class), abx_class != ""
  ) %>%
  semi_join(cohort_primary %>% distinct(infant_study_id), by = "infant_study_id") %>%
  distinct(abx_class) %>%
  pull(abx_class) %>%
  as.character()

stopifnot(
  "Reviewer 3 target antibiotic class label(s) not found in abx_class" =
    all(TARGET_CLASSES_R3 %in% class_levels_r3)
)

build_r3_class_events <- function(class_names) {
  analytic_events %>%
    filter(
      analyzable_event_primary == TRUE,
      abx_class %in% class_names,
      !is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX
    ) %>%
    semi_join(cohort_primary %>% distinct(infant_study_id), by = "infant_study_id") %>%
    distinct(infant_study_id, abx_days)
}

COMPOSITE_CLASS_R3 <- "Aminoglycoside or vancomycin"

TARGET_OUTCOMES_R3 <- c(COMPOSITE_CLASS_R3, TARGET_CLASSES_R3)

target_class_events_r3 <- c(
  set_names(list(build_r3_class_events(TARGET_CLASSES_R3)), COMPOSITE_CLASS_R3),
  set_names(map(TARGET_CLASSES_R3, build_r3_class_events), TARGET_CLASSES_R3)
)
print(imap_dfr(target_class_events_r3, ~ tibble(
  antibiotic_class = .y,
  n_event_days = nrow(.x),
  n_infants_with_event = n_distinct(.x$infant_study_id)
)), n = Inf)

add_calendar_3month_r3 <- function(df) {
  df %>%
    mutate(
      calendar_3month = factor(paste0(
        format(calendar_date_start, "%Y"), "_Q",
        (as.integer(format(calendar_date_start, "%m")) - 1L) %/% 3L + 1L
      ))
    )
}

COV_CLASS_ASSOC_R3 <- paste(
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED +",
  "tot_covid + gbs_prior + utd_status + lb_prior"
)

formula_class_assoc_r3 <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + flu_prior + tdap_prior +",
  "strata(calendar_3month) +", COV_CLASS_ASSOC_R3
))

fit_class_assoc_r3 <- function(class_label,
                               period = c("full_year", "rsv_season")) {
  period <- match.arg(period)

  cp <- build_cp(target_class_events_r3[[class_label]], blocks = blocks_primary)

  md <- switch(
    period,
    full_year  = clean_and_restrict(cp, restrict_to_season = FALSE),
    rsv_season = clean_and_restrict(cp, restrict_to_season = TRUE)
  ) %>%
    add_calendar_3month_r3()

  n_events <- sum(md$event)
  n_infants <- n_distinct(md$infant_study_id)

  if (n_events < MIN_EVENTS_TO_FIT) {
    return(tibble(
      antibiotic_class = class_label,
      period = period,
      term = c("flu_prior", "tdap_prior"),
      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p.value = NA_real_,
      n_infants = n_infants,
      n_events = n_events,
      fit_status = paste0("Skipped: only ", n_events, " events")
    ))
  }

  fit <- coxph(
    formula_class_assoc_r3,
    data = md,
    cluster = mom_study_id,
    ties = "efron",
    model = TRUE
  )

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("flu_prior", "tdap_prior")) %>%
    transmute(
      antibiotic_class = class_label,
      period = period,
      term,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p.value,
      n_infants = n_infants,
      n_events = n_events,
      fit_status = "OK"
    )
}

r3_flu_tdap_class_assoc <- crossing(
  antibiotic_class = TARGET_OUTCOMES_R3,
  period = c("full_year", "rsv_season")
) %>%
  pmap_dfr(~ fit_class_assoc_r3(..1, ..2)) %>%
  mutate(
    maternal_vaccine = recode(
      term,
      flu_prior  = "Maternal influenza vaccination",
      tdap_prior = "Maternal Tdap vaccination"
    ),
    HR_95CI = if_else(
      is.na(HR), NA_character_,
      sprintf("%.2f (%.2f to %.2f)", HR, CI_low, CI_high)
    ),
    period_label = recode(
      period,
      full_year = "Full year, days 1-180",
      rsv_season = "RSV season, days 1-180"
    )
  ) %>%
  mutate(antibiotic_class = factor(antibiotic_class, levels = TARGET_OUTCOMES_R3)) %>%
  arrange(antibiotic_class, factor(period, levels = c("full_year", "rsv_season"))) %>%
  mutate(antibiotic_class = as.character(antibiotic_class)) %>%
  select(
    antibiotic_class, period, period_label, maternal_vaccine,
    n_infants, n_events, HR, CI_low, CI_high, HR_95CI, p.value, fit_status
  )
print(r3_flu_tdap_class_assoc, n = Inf, width = Inf)

write_csv(r3_ve_results, "reviewer3_VE_sensitivity_long.csv")
write_csv(r3_ve_table_wide, "reviewer3_VE_sensitivity_table.csv")
write_csv(r3_ve_change, "reviewer3_VE_change_vs_primary.csv")
write_csv(r3_flu_tdap_class_assoc, "reviewer3_flu_tdap_antibiotic_class_associations.csv")
