

library(gtsummary)
library(gt)

stopifnot(exists("analytic_infants"), exists("CANONICAL_N"))

table1_cohort <- analytic_infants %>%
  filter(has_complete_covariates)
print(CANONICAL_N$complete_case %>% select(exposure, n_infants_cohort, total, pct))

stopifnot(
  nrow(table1_cohort) == sum(CANONICAL_N$complete_case$n_infants_cohort),
  !any(is.na(table1_cohort$exposure_group)),

  nrow(table1_cohort) == n_distinct(model_data_fullyear$all_cause$infant_study_id),
  nrow(table1_cohort) >= n_distinct(model_data_primary$all_cause$infant_study_id)
)

UTD_LANDMARK_DAY <- 90

table1_cohort <- table1_cohort %>%
  build_utd_at(UTD_LANDMARK_DAY) %>%
  select(-utd_hepb, -utd_rota, -utd_dtap, -utd_hib, -utd_pcv, -utd_ipv)

strip_leading_code <- function(x) stringr::str_remove(as.character(x), "^\\d+:\\s*")

diagnose_binary <- function(x, label) {
  print(table(x, useNA = "ifany"))
}

diagnose_binary(table1_cohort$serious_newborn_comorbidity, "serious_newborn_comorbidity")
diagnose_binary(table1_cohort$gbs_prior,   "gbs_prior")
diagnose_binary(table1_cohort$flu_prior,   "flu_prior")
diagnose_binary(table1_cohort$tdap_prior,  "tdap_prior")

as_yesno <- function(x) {
  if (is.logical(x)) {
    lvl <- dplyr::if_else(x, "Yes", "No")
  } else if (is.numeric(x)) {
    lvl <- dplyr::case_when(x == 1 ~ "Yes", x == 0 ~ "No", TRUE ~ NA_character_)
  } else {
    x_chr <- toupper(trimws(as.character(x)))
    lvl <- dplyr::case_when(
      x_chr %in% c("YES", "Y", "1", "TRUE")  ~ "Yes",
      x_chr %in% c("NO", "N", "0", "FALSE")  ~ "No",
      TRUE ~ NA_character_
    )
  }
  factor(lvl, levels = c("No", "Yes"))
}

table1_cohort <- table1_cohort %>%
  mutate(
    GENDER = forcats::fct_drop(factor(dplyr::if_else(GENDER == "" | is.na(GENDER),
                                                     "Unknown", as.character(GENDER)))),
    race_cat_infant = forcats::fct_drop(forcats::fct_relevel(
      factor(race_cat_infant),
      "Non-Hispanic Asian/Pacific Islander", "Hispanic",
      "Non-Hispanic Black", "Non-Hispanic White", "Other/Unknown"
    )),
    race_cat_mom = forcats::fct_drop(forcats::fct_relevel(
      factor(race_cat_mom),
      "Non-Hispanic Asian/Pacific Islander", "Hispanic",
      "Non-Hispanic Black", "Non-Hispanic White", "Other/Unknown"
    )),
    smoking = forcats::fct_drop(forcats::fct_relevel(
      factor(dplyr::if_else(is.na(smoking), "Unknown", as.character(smoking))),
      "Current smoker", "Former smoker", "Never smoker", "Unknown"
    )),
    median_income_cat = forcats::fct_drop(forcats::fct_relevel(
      factor(strip_leading_code(median_income_cat)),
      "< $40,000", "$40,000 - $79,999", ">= $80,000", "Unknown"
    )),
    utd_status = forcats::fct_drop(utd_status),
    serious_newborn_comorbidity = as_yesno(serious_newborn_comorbidity),
    gbs_prior  = as_yesno(gbs_prior),
    flu_prior  = as_yesno(flu_prior),
    tdap_prior = as_yesno(tdap_prior)
  )

table1_vars <- c(
  "gestage_imp", "apgar_1min", "GENDER", "utd_status", "race_cat_infant",
  "serious_newborn_comorbidity", "age_index", "race_cat_mom", "bmi_prior",
  "smoking", "median_income_cat", "charlson_Wt", "gbs_prior",
  "flu_prior", "tdap_prior", "OUTPT", "INPT", "ED", "tot_covid"
)

TABLE1_COLUMN_ORDER <- c("Maternal RSV vaccine only", "Nirsevimab only", "Neither")

table1_gt <- table1_cohort %>%
  mutate(exposure_group = factor(
    unname(EXPOSURE_GROUP_DISPLAY[as.character(exposure_group)]),
    levels = unname(EXPOSURE_GROUP_DISPLAY[TABLE1_COLUMN_ORDER])
  )) %>%
  select(exposure_group, all_of(table1_vars)) %>%
  tbl_summary(
    by = exposure_group,
    type = list(
      charlson_Wt ~ "continuous", tot_covid ~ "continuous",
      serious_newborn_comorbidity ~ "dichotomous",
      gbs_prior ~ "dichotomous", flu_prior ~ "dichotomous", tdap_prior ~ "dichotomous"
    ),
    value = list(
      serious_newborn_comorbidity ~ "Yes",
      gbs_prior ~ "Yes", flu_prior ~ "Yes", tdap_prior ~ "Yes"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),

    missing = "ifany",
    missing_text = "Unknown, n (%)",
    missing_stat = "{N_miss} ({p_miss}%)",
    label = list(
      gestage_imp                  ~ "Gestational age at birth (days)",
      apgar_1min                   ~ "1-minute Apgar score",
      GENDER                       ~ "Infant sex",
      utd_status                   ~ paste0("Up to date on routine vaccines at day ", UTD_LANDMARK_DAY),
      race_cat_infant              ~ "Infant race/ethnicity",
      serious_newborn_comorbidity  ~ "Serious newborn comorbidity",
      age_index                    ~ "Maternal age at delivery (years)",
      race_cat_mom                 ~ "Maternal race/ethnicity",
      bmi_prior                    ~ "Maternal pre-pregnancy BMI",
      smoking                      ~ "Maternal smoking status",
      median_income_cat            ~ "Median household income category",
      charlson_Wt                  ~ "Maternal weighted Charlson comorbidity index",
      gbs_prior                    ~ "Pre-delivery GBS screening",
      flu_prior                    ~ "Maternal influenza vaccine, 1yr pre-delivery",
      tdap_prior                   ~ "Maternal Tdap vaccine, 1yr pre-delivery",
      OUTPT                        ~ "Maternal outpatient visit count, 1yr pre-delivery",
      INPT                         ~ "Maternal inpatient visit count, 1yr pre-delivery",
      ED                           ~ "Maternal ED visit count, 1yr pre-delivery",
      tot_covid                    ~ "Total maternal COVID-19 vaccine doses"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  add_n() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_caption(paste0(
    "**Table 1. Baseline characteristics of the never-both cohort with ",
    "complete covariate data, by exposure group**"
  )) %>%
  bold_labels()

print(table1_gt)

table1_gt %>%
  as_gt() %>%
  gt::tab_source_note(source_note = gt::md(paste0(
    "Exposure groups are mutually exclusive and are assigned at the infant ",
    "level; infants who received both maternal RSV vaccine and nirsevimab ",
    "are excluded from this cohort by design. Counts therefore sum to the ",
    "column total, and differ from the exposure-specific infant counts in ",
    "Table S1 and in the forest plots, which are model-specific and reflect ",
    "time-varying nirsevimab exposure. Up-to-date vaccination status is ",
    "evaluated at day ", UTD_LANDMARK_DAY, " of follow-up."
  ))) %>%
  gt::tab_source_note(source_note = gt::md(paste0(
    "Population: never-both cohort infants with non-missing values on every ",
    "adjustment covariate (n = ", scales::comma(nrow(table1_cohort)), "). ",
    "This population is defined independently of RSV season restriction and ",
    "does not change with the RSV_SEASONS window; it therefore differs from, ",
    "and is larger than or equal to, the season-restricted N reported for the ",
    "primary Cox model."
  ))) %>%
  gt::gtsave("table1_baseline_characteristics2.docx")
