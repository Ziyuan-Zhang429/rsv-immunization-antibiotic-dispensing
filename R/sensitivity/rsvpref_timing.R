

library(dplyr)
library(gtsummary)
library(gt)
library(forcats)
library(stringr)

if (!exists("eligible_cohort")) {
  stop(paste0(
    "Object `eligible_cohort` does not exist. From a fresh R session, first run ",
    "AG_Cox_postindex_antibiotic_calendar_month_strata.R from the beginning ",
    "through the end of Phase 3.2 (sibling-based reclassification), then run this file."
  ))
}

required_cols <- c(
  "rsv_history", "rsv_valid", "sibling_match_confirmed",
  "gestage_at_vax", "days_vax_to_birth",
  "bronchopulm_dx", "immunodef_dx", "downsyn_dx",
  "gestage_imp", "apgar_1min", "GENDER", "race_cat_infant",
  "age_index", "race_cat_mom", "bmi_prior", "smoking",
  "median_income_cat", "charlson_Wt", "gbs_prior", "flu_prior",
  "tdap_prior", "OUTPT", "INPT", "ED", "tot_covid"
)
missing_cols <- setdiff(required_cols, names(eligible_cohort))
if (length(missing_cols) > 0) {
  stop("Missing required columns in eligible_cohort: ", paste(missing_cols, collapse = ", "))
}

if (!exists("GESTAGE_AT_VAX_MIN")) GESTAGE_AT_VAX_MIN <- 224
if (!exists("GESTAGE_AT_VAX_MAX")) GESTAGE_AT_VAX_MAX <- 252
if (!exists("MIN_DAYS_VAX_TO_BIRTH")) MIN_DAYS_VAX_TO_BIRTH <- 14

rsv_timing_cohort <- eligible_cohort %>%
  filter(rsv_history, !sibling_match_confirmed) %>%
  mutate(
    rsv_timing_group = factor(
      if_else(rsv_valid, "Eligible RSVpreF timing", "Ineligible RSVpreF timing"),
      levels = c("Eligible RSVpreF timing", "Ineligible RSVpreF timing")
    ),
    serious_newborn_comorbidity = as.integer(
      bronchopulm_dx == 1 | immunodef_dx == 1 | downsyn_dx == 1
    ),
    eligible_gestage = !is.na(gestage_at_vax) &
      gestage_at_vax >= GESTAGE_AT_VAX_MIN &
      gestage_at_vax <= GESTAGE_AT_VAX_MAX,
    eligible_14day = !is.na(days_vax_to_birth) &
      days_vax_to_birth >= MIN_DAYS_VAX_TO_BIRTH,
    ineligible_reason = case_when(
      rsv_valid ~ "Eligible",
      !eligible_gestage & !eligible_14day ~ "Outside 32-36 weeks and <14 days before delivery",
      !eligible_gestage ~ "Outside 32-36 weeks",
      !eligible_14day ~ "<14 days before delivery",
      TRUE ~ "Other/unknown"
    )
  )
print(rsv_timing_cohort %>%
        count(rsv_timing_group) %>%
        mutate(pct = 100 * n / sum(n)))
print(rsv_timing_cohort %>%
        filter(rsv_timing_group == "Ineligible RSVpreF timing") %>%
        count(ineligible_reason) %>%
        mutate(pct = 100 * n / sum(n)))

write.csv(
  rsv_timing_cohort %>% count(rsv_timing_group) %>% mutate(pct = 100 * n / sum(n)),
  "selection_bias_B_counts.csv", row.names = FALSE
)
write.csv(
  rsv_timing_cohort %>%
    filter(rsv_timing_group == "Ineligible RSVpreF timing") %>%
    count(ineligible_reason) %>% mutate(pct = 100 * n / sum(n)),
  "selection_bias_B_ineligible_reasons.csv", row.names = FALSE
)

strip_leading_code <- function(x) {
  stringr::str_remove(as.character(x), "^\\d+:\\s*")
}

as_yesno <- function(x) {
  if (is.logical(x)) {
    lvl <- dplyr::if_else(x, "Yes", "No", missing = NA_character_)
  } else if (is.numeric(x)) {
    lvl <- dplyr::case_when(x == 1 ~ "Yes", x == 0 ~ "No", TRUE ~ NA_character_)
  } else {
    x_chr <- toupper(trimws(as.character(x)))
    lvl <- dplyr::case_when(
      x_chr %in% c("YES", "Y", "1", "TRUE") ~ "Yes",
      x_chr %in% c("NO", "N", "0", "FALSE") ~ "No",
      TRUE ~ NA_character_
    )
  }
  factor(lvl, levels = c("No", "Yes"))
}

rsv_timing_cohort <- rsv_timing_cohort %>%
  mutate(

    GENDER = forcats::fct_drop(
      factor(
        if_else(
          GENDER == "" | is.na(GENDER),
          "Unknown",
          as.character(GENDER)
        )
      )
    ),

    race_cat_infant = forcats::fct_drop(
      forcats::fct_relevel(
        factor(race_cat_infant),
        "Asian/PI",
        "Hispanic",
        "Non-hispanic Black",
        "Non-hispanic White",
        "Other/Unknown"
      )
    ),

    race_cat_mom = forcats::fct_drop(
      forcats::fct_relevel(
        factor(race_cat_mom),
        "Asian/PI",
        "Hispanic",
        "Non-hispanic Black",
        "Non-hispanic White",
        "Other/Unknown"
      )
    ),

    smoking = forcats::fct_drop(
      forcats::fct_relevel(
        factor(
          if_else(
            is.na(smoking),
            "Unknown",
            as.character(smoking)
          )
        ),
        "Current smoker",
        "Former smoker",
        "Never smoker",
        "Unknown"
      )
    ),

    median_income_cat = forcats::fct_drop(
      forcats::fct_relevel(
        factor(strip_leading_code(median_income_cat)),
        "< $40,000",
        "$40,000 - $79,999",
        ">= $80,000"
      )
    ),

    serious_newborn_comorbidity =
      as_yesno(serious_newborn_comorbidity),

    gbs_prior  = as_yesno(gbs_prior),
    flu_prior  = as_yesno(flu_prior),
    tdap_prior = as_yesno(tdap_prior)
  )

selection_vars <- c(
  "gestage_imp", "apgar_1min", "GENDER", "race_cat_infant",
  "serious_newborn_comorbidity", "age_index", "race_cat_mom",
  "bmi_prior", "smoking", "median_income_cat", "charlson_Wt",
  "gbs_prior", "flu_prior", "tdap_prior", "OUTPT", "INPT", "ED", "tot_covid"
)

table_selection_rsv <- rsv_timing_cohort %>%
  select(rsv_timing_group, all_of(selection_vars)) %>%
  tbl_summary(
    by = rsv_timing_group,
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
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "ifany",
    missing_text = "Unknown, n (%)",
    missing_stat = "{N_miss} ({p_miss}%)",
    label = list(
      gestage_imp ~ "Gestational age at birth (days)",
      apgar_1min ~ "1-minute Apgar score",
      GENDER ~ "Infant sex",
      race_cat_infant ~ "Infant race/ethnicity",
      serious_newborn_comorbidity ~ "Serious newborn comorbidity",
      age_index ~ "Maternal age at delivery (years)",
      race_cat_mom ~ "Maternal race/ethnicity",
      bmi_prior ~ "Maternal pre-pregnancy BMI",
      smoking ~ "Maternal smoking status",
      median_income_cat ~ "Median household income category",
      charlson_Wt ~ "Maternal weighted Charlson comorbidity index",
      gbs_prior ~ "Pre-delivery GBS screening",
      flu_prior ~ "Maternal influenza vaccine, 1yr pre-delivery",
      tdap_prior ~ "Maternal Tdap vaccine, 1yr pre-delivery",
      OUTPT ~ "Maternal outpatient visit count, 1yr pre-delivery",
      INPT ~ "Maternal inpatient visit count, 1yr pre-delivery",
      ED ~ "Maternal ED visit count, 1yr pre-delivery",
      tot_covid ~ "Total maternal COVID-19 vaccine doses"
    )
  ) %>%
  add_difference(test = everything() ~ "smd") %>%
  modify_header(label ~ "**Characteristic**", estimate ~ "**SMD**") %>%
  modify_caption("**Table SXX. Characteristics of maternal RSVpreF recipients by eligibility of vaccination timing**") %>%
  bold_labels()

print(table_selection_rsv)

table_selection_rsv %>%
  as_gt() %>%
  gt::tab_source_note(source_note = gt::md(
    paste0(
      "Eligible maternal RSVpreF timing was defined as vaccination at 32-36 weeks' gestation ",
      "and at least 14 days before delivery. Ineligible timing included maternal RSVpreF recipients ",
      "who did not meet one or both timing criteria. Records identified as belonging to a later ",
      "sibling's pregnancy were excluded from this comparison."
    )
  )) %>%
  gt::tab_source_note(source_note = gt::md(
    "Values are mean (SD) for continuous variables and n (%) for categorical variables. SMD = standardized mean difference."
  )) %>%
  gt::gtsave("Table_SXX_selection_bias_RSVpreF_timing.docx")
