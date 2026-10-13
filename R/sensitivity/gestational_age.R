

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(survival)
  library(broom)
  library(readr)
  library(tibble)
})

BASE_SCRIPT <- "R/primary/cox_primary.R"

if (!file.exists(BASE_SCRIPT)) {
  stop("Cannot find ", BASE_SCRIPT, " in the working directory.")
}
source(BASE_SCRIPT)

stopifnot(
  exists("cp_primary"),
  exists("clean_and_restrict"),
  exists("fit_ag"),
  exists("tidy_ve"),
  exists("COVARIATE_TERMS"),
  exists("primary_ve_table"),
  exists("sens_ve_table_fullyear"),
  exists("out_of_season_ve_table"),
  exists("model_data_primary"),
  exists("model_data_fullyear"),
  exists("model_data_out_of_season")
)

if (!"gestage_imp" %in% names(cp_primary$all_cause)) {
  stop("gestage_imp is not present in cp_primary$all_cause.")
}

original_covariate_terms <- COVARIATE_TERMS
original_primary <- primary_ve_table
original_fullyear <- sens_ve_table_fullyear
original_oos <- out_of_season_ve_table
add_preterm <- function(dat) {
  dat %>%
    mutate(
      preterm_birth = case_when(
        is.na(gestage_imp) ~ NA_integer_,
        gestage_imp < 259  ~ 1L,
        TRUE               ~ 0L
      )
    )
}

cp_ga <- map(cp_primary, add_preterm)

ga_infant_check <- cp_ga$all_cause %>%
  distinct(infant_study_id, gestage_imp, preterm_birth)

ga_distribution <- ga_infant_check %>%
  summarise(
    n_infants = n(),
    n_missing_ga = sum(is.na(gestage_imp)),
    n_preterm_32_to_lt37 = sum(preterm_birth == 1L, na.rm = TRUE),
    pct_preterm_32_to_lt37 = 100 * mean(preterm_birth == 1L, na.rm = TRUE),
    n_term_ge37 = sum(preterm_birth == 0L, na.rm = TRUE),
    pct_term_ge37 = 100 * mean(preterm_birth == 0L, na.rm = TRUE),
    min_ga_days = min(gestage_imp, na.rm = TRUE),
    median_ga_days = median(gestage_imp, na.rm = TRUE),
    max_ga_days = max(gestage_imp, na.rm = TRUE)
  )
print(ga_distribution)
COVARIATE_TERMS_GA <- paste(
  "preterm_birth +",
  original_covariate_terms
)

COVARIATE_TERMS <- COVARIATE_TERMS_GA

ag_formula_ga <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_month) +",
  COVARIATE_TERMS_GA
))

model_data_primary_ga <- map(
  cp_ga,
  ~ clean_and_restrict(.x, restrict_to_season = TRUE)
)

model_data_fullyear_ga <- map(
  cp_ga,
  ~ clean_and_restrict(.x, restrict_to_season = FALSE)
)

model_data_out_of_season_ga <- map(
  cp_ga,
  ~ clean_and_restrict(.x, keep_out_of_season_only = TRUE)
)

sample_size_comparison <- bind_rows(
  tibble(
    period = "RSV seasons",
    original_n = n_distinct(model_data_primary$all_cause$infant_study_id),
    preterm_adjusted_n = n_distinct(model_data_primary_ga$all_cause$infant_study_id)
  ),
  tibble(
    period = "Full year",
    original_n = n_distinct(model_data_fullyear$all_cause$infant_study_id),
    preterm_adjusted_n = n_distinct(model_data_fullyear_ga$all_cause$infant_study_id)
  ),
  tibble(
    period = "Non-RSV seasons",
    original_n = n_distinct(model_data_out_of_season$all_cause$infant_study_id),
    preterm_adjusted_n = n_distinct(model_data_out_of_season_ga$all_cause$infant_study_id)
  )
) %>%
  mutate(n_difference = preterm_adjusted_n - original_n)
print(sample_size_comparison)

cox_primary_ga <- map(
  model_data_primary_ga,
  ~ fit_ag(.x, formula = ag_formula_ga)
)
primary_ve_table_ga <- imap_dfr(
  cox_primary_ga,
  ~ tidy_ve(.x, .y)
)

cox_fullyear_ga <- map(
  model_data_fullyear_ga,
  ~ fit_ag(.x, formula = ag_formula_ga)
)
sens_ve_table_fullyear_ga <- imap_dfr(
  cox_fullyear_ga,
  ~ tidy_ve(.x, .y)
)

cox_out_of_season_ga <- map(
  model_data_out_of_season_ga,
  ~ fit_ag(.x, formula = ag_formula_ga)
)
out_of_season_ve_table_ga <- imap_dfr(
  cox_out_of_season_ga,
  ~ tidy_ve(.x, .y)
)

make_comparison <- function(original, ga_adjusted, period_label) {
  original %>%
    select(
      outcome, exposure,
      original_HR = HR,
      original_CI_low = CI_low,
      original_CI_high = CI_high,
      original_VE_pct = VE_pct,
      original_VE_CI_low = VE_CI_low,
      original_VE_CI_high = VE_CI_high
    ) %>%
    left_join(
      ga_adjusted %>%
        select(
          outcome, exposure,
          preterm_adjusted_HR = HR,
          preterm_adjusted_CI_low = CI_low,
          preterm_adjusted_CI_high = CI_high,
          preterm_adjusted_VE_pct = VE_pct,
          preterm_adjusted_VE_CI_low = VE_CI_low,
          preterm_adjusted_VE_CI_high = VE_CI_high
        ),
      by = c("outcome", "exposure")
    ) %>%
    mutate(
      period = period_label,
      change_in_VE_percentage_points =
        preterm_adjusted_VE_pct - original_VE_pct,
      .before = 1
    )
}

ga_adjustment_comparison_all <- bind_rows(
  make_comparison(original_primary, primary_ve_table_ga, "RSV seasons"),
  make_comparison(original_fullyear, sens_ve_table_fullyear_ga, "Full year"),
  make_comparison(original_oos, out_of_season_ve_table_ga, "Non-RSV seasons")
)

ga_adjustment_comparison_all_cause <- ga_adjustment_comparison_all %>%
  filter(
    outcome == "all_cause",
    exposure %in% c(
      "exposure_arma_maternal_vax_only",
      "exposure_armb_nirsevimab_only"
    )
  )
print(ga_adjustment_comparison_all_cause, n = Inf, width = Inf)
COVARIATE_TERMS <- original_covariate_terms

md_all <- clean_and_restrict(
  cp_ga$all_cause,
  restrict_to_season = TRUE
) %>%
  mutate(
    ga_group = factor(
      preterm_birth,
      levels = c(1, 0),
      labels = c("Preterm (32 to <37 weeks)", "Term (>=37 weeks)")
    )
  )

if (any(is.na(md_all$preterm_birth))) {
  warning("Rows with missing gestational age are excluded from the subgroup analysis.")
  md_all <- md_all %>% filter(!is.na(preterm_birth))
}

ga_group_counts <- md_all %>%
  distinct(infant_study_id, ga_group, preterm_birth) %>%
  count(ga_group, name = "n_infants") %>%
  mutate(percent = 100 * n_infants / sum(n_infants))

ga_group_events <- md_all %>%
  group_by(ga_group) %>%
  summarise(
    n_infants = n_distinct(infant_study_id),
    antibiotic_courses = sum(event, na.rm = TRUE),
    person_days = sum(tstop - tstart, na.rm = TRUE),
    .groups = "drop"
  )
print(ga_group_counts)
print(ga_group_events)

strat_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_month) +",
  original_covariate_terms
))

fit_one_ga_group <- function(dat, label) {
  fit <- coxph(
    strat_formula,
    data = dat,
    cluster = mom_study_id,
    ties = "efron",
    model = TRUE
  )

  tidy_ve(fit, "all_cause") %>%
    mutate(
      ga_group = label,
      n_infants = n_distinct(dat$infant_study_id),
      antibiotic_courses = sum(dat$event, na.rm = TRUE)
    ) %>%
    select(
      ga_group, outcome, exposure,
      HR, CI_low, CI_high,
      VE_pct, VE_CI_low, VE_CI_high, p.value,
      n_infants, antibiotic_courses
    )
}

stratified_ve <- bind_rows(
  fit_one_ga_group(
    md_all %>% filter(preterm_birth == 1L),
    "Preterm (32 to <37 weeks)"
  ),
  fit_one_ga_group(
    md_all %>% filter(preterm_birth == 0L),
    "Term (>=37 weeks)"
  )
)
print(stratified_ve, n = Inf, width = Inf)

interaction_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm * preterm_birth +",
  "strata(calendar_month) +",
  original_covariate_terms
))

interaction_fit <- coxph(
  interaction_formula,
  data = md_all,
  cluster = mom_study_id,
  ties = "efron",
  model = TRUE
)

interaction_tidy <- broom::tidy(
  interaction_fit,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(grepl("^exposure_arm.*:preterm_birth$", term)) %>%
  transmute(
    interaction_term = term,
    ratio_of_HRs = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p_value_interaction = p.value
  )

int_names <- names(coef(interaction_fit))[
  grepl("^exposure_arm.*:preterm_birth$", names(coef(interaction_fit)))
]

if (length(int_names) > 0) {
  b <- coef(interaction_fit)[int_names]
  V <- vcov(interaction_fit)[int_names, int_names, drop = FALSE]

  joint_stat <- tryCatch(
    as.numeric(t(b) %*% solve(V, b)),
    error = function(e) NA_real_
  )

  joint_interaction_test <- tibble(
    test = "Joint Wald test: exposure_arm x preterm_birth",
    df = length(int_names),
    chi_square = joint_stat,
    p_value = ifelse(
      is.na(joint_stat),
      NA_real_,
      pchisq(joint_stat, df = length(int_names), lower.tail = FALSE)
    )
  )
} else {
  joint_interaction_test <- tibble(
    test = "Joint Wald test: exposure_arm x preterm_birth",
    df = 0L,
    chi_square = NA_real_,
    p_value = NA_real_
  )
}
print(interaction_tidy, n = Inf, width = Inf)
print(joint_interaction_test)

ga_exposure_diagnostics <- md_all %>%
  group_by(ga_group, exposure_arm) %>%
  summarise(
    n_infants_contributing_person_time = n_distinct(infant_study_id),
    antibiotic_courses = sum(event, na.rm = TRUE),
    person_days = sum(tstop - tstart, na.rm = TRUE),
    .groups = "drop"
  )
print(ga_exposure_diagnostics, n = Inf)
fit_sanity <- coxph(
  strat_formula,
  data = md_all,
  cluster = mom_study_id,
  ties = "efron",
  model = TRUE
)

sanity_ve <- tidy_ve(fit_sanity, "all_cause")

original_primary_allcause <- original_primary %>%
  filter(outcome == "all_cause")

sanity_comparison <- sanity_ve %>%
  select(
    exposure,
    sanity_HR = HR,
    sanity_VE = VE_pct
  ) %>%
  left_join(
    original_primary_allcause %>%
      select(
        exposure,
        original_HR = HR,
        original_VE = VE_pct
      ),
    by = "exposure"
  ) %>%
  mutate(
    HR_difference = sanity_HR - original_HR,
    VE_difference = sanity_VE - original_VE
  )
print(sanity_comparison, n = Inf, width = Inf)

write_csv(
  ga_distribution,
  "reviewer4_GA_distribution.csv"
)

write_csv(
  sample_size_comparison,
  "reviewer4_GA_adjustment_sample_sizes.csv"
)

write_csv(
  ga_adjustment_comparison_all_cause,
  "reviewer4_GA_adjustment_all_cause.csv"
)

write_csv(
  ga_adjustment_comparison_all,
  "reviewer4_GA_adjustment_all_outcomes.csv"
)

write_csv(
  stratified_ve,
  "reviewer4_GA_preterm_term_stratified_VE.csv"
)

write_csv(
  ga_group_events,
  "reviewer4_GA_preterm_term_sample_events.csv"
)

write_csv(
  interaction_tidy,
  "reviewer4_GA_preterm_term_interaction.csv"
)

write_csv(
  joint_interaction_test,
  "reviewer4_GA_preterm_term_joint_interaction.csv"
)

write_csv(
  ga_exposure_diagnostics,
  "reviewer4_GA_preterm_term_exposure_diagnostics.csv"
)

write_csv(
  sanity_comparison,
  "reviewer4_GA_sanity_check.csv"
)

COVARIATE_TERMS <- original_covariate_terms
