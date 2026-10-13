

library(dplyr)
library(tidyr)
library(gt)
library(scales)

stopifnot(
  exists("cp_primary"),
  "all_cause" %in% names(cp_primary),
  exists("COVARIATE_TERMS"),
  exists("analytic_infants")
)

cp0 <- cp_primary$all_cause %>%
  filter(!is.na(exposure_arm))

eligible_ids <- analytic_infants %>%
  distinct(infant_study_id)

cp0 <- cp0 %>%
  semi_join(eligible_ids, by = "infant_study_id")

stopifnot(n_distinct(cp0$infant_study_id) == 47868)

model_predictors <- all.vars(as.formula(paste("~", COVARIATE_TERMS)))
print(model_predictors)

stopifnot(all(model_predictors %in% names(cp0)))

pred_missing_long <- cp0 %>%
  transmute(
    infant_study_id,
    across(all_of(model_predictors), is.na, .names = "miss__{.col}")
  ) %>%
  group_by(infant_study_id) %>%
  summarise(
    across(starts_with("miss__"), any),
    .groups = "drop"
  )

unknown_flags <- cp0 %>%
  group_by(infant_study_id) %>%
  summarise(
    unk__GENDER = any(GENDER == "", na.rm = TRUE),
    unk__median_income_cat = any(median_income_cat == "4: Unknown", na.rm = TRUE),
    unk__smoking = any(smoking == "Unknown", na.rm = TRUE),
    unk__bmi_cat = any(bmi_cat == "5: Unknown", na.rm = TRUE),
    .groups = "drop"
  )

flags <- eligible_ids %>%
  left_join(pred_missing_long, by = "infant_study_id") %>%
  left_join(unknown_flags, by = "infant_study_id") %>%
  mutate(
    across(starts_with("miss__"), ~ replace_na(.x, FALSE)),
    across(starts_with("unk__"),  ~ replace_na(.x, FALSE))
  )

flag_for <- function(var) {
  nm <- paste0("miss__", var)
  if (!nm %in% names(flags)) return(rep(FALSE, nrow(flags)))
  flags[[nm]]
}

reason_flags <- tibble(
  infant_study_id = flags$infant_study_id,

  `Infant sex` =
    flag_for("GENDER") | flags$unk__GENDER,

  `Infant race/ethnicity` =
    flag_for("race_cat_infant"),

  `Serious newborn comorbidity` =
    flag_for("serious_newborn_comorbidity"),

  `1-minute Apgar score` =
    flag_for("apgar_1min"),

  `Maternal age at delivery` =
    flag_for("age_index"),

  `Maternal weighted Charlson comorbidity index` =
    flag_for("charlson_Wt"),

  `Maternal body mass index (pre-pregnancy)` =
    flag_for("bmi_prior") | flags$unk__bmi_cat,

  `Maternal cigarette smoking` =
    flag_for("smoking") | flags$unk__smoking,

  `Neighborhood median household income` =
    flag_for("median_income_cat") | flags$unk__median_income_cat,

  `Maternal outpatient visit count, 1yr pre-delivery` =
    flag_for("OUTPT"),

  `Maternal inpatient admission count, 1yr pre-delivery` =
    flag_for("INPT"),

  `Maternal ED visit count, 1yr pre-delivery` =
    flag_for("ED"),

  `Maternal influenza vaccine, 1yr pre-delivery` =
    flag_for("flu_prior"),

  `Maternal Tdap vaccine, 1yr pre-delivery` =
    flag_for("tdap_prior"),

  `Total maternal COVID-19 vaccine doses` =
    flag_for("tot_covid"),

  `Pre-delivery group B Streptococcus screening` =
    flag_for("gbs_prior"),

  `Age-appropriate routine vaccination status` =
    flag_for("utd_status"),

  `Prior live birth` =
    flag_for("lb_prior")
)

reason_cols <- setdiff(names(reason_flags), "infant_study_id")

reason_flags <- reason_flags %>%
  mutate(any_missing_or_unknown = if_any(all_of(reason_cols), identity))

validation <- analytic_infants %>%
  select(infant_study_id, has_complete_covariates) %>%
  left_join(
    reason_flags %>% select(infant_study_id, any_missing_or_unknown),
    by = "infant_study_id"
  ) %>%
  mutate(
    any_missing_or_unknown = replace_na(any_missing_or_unknown, FALSE),
    expected_excluded = !has_complete_covariates
  )
print(table(
  reconstructed_excluded = validation$any_missing_or_unknown,
  canonical_excluded = validation$expected_excluded
))

n_canonical <- sum(validation$expected_excluded)
n_reconstructed <- sum(validation$any_missing_or_unknown)
n_discordant <- sum(
  validation$expected_excluded != validation$any_missing_or_unknown
)
stopifnot(n_canonical == 218)

if (n_discordant > 0) {
  warning(
    paste0(
      "The reconstructed variable-level rules do not yet exactly reproduce ",
      "has_complete_covariates (", n_discordant, " discordant infants). ",
      "Inspect `discordant_infants` before using the table."
    )
  )
}

discordant_infants <- validation %>%
  filter(expected_excluded != any_missing_or_unknown)

missing_table_data <- reason_flags %>%
  summarise(across(all_of(reason_cols), sum)) %>%
  pivot_longer(
    everything(),
    names_to = "characteristic",
    values_to = "n_missing"
  ) %>%
  mutate(
    pct_missing = 100 * n_missing / nrow(analytic_infants),
    `Missing/unknown, n (%)` = sprintf(
      "%s (%.2f%%)",
      comma(n_missing),
      pct_missing
    )
  )

missing_table_nonzero <- missing_table_data %>%
  filter(n_missing > 0)

any_row <- tibble(
  characteristic = "Any required covariate",
  n_missing = n_canonical,
  pct_missing = 100 * n_canonical / nrow(analytic_infants),
  `Missing/unknown, n (%)` = sprintf(
    "%s (%.2f%%)",
    comma(n_canonical),
    100 * n_canonical / nrow(analytic_infants)
  )
)

missing_table_with_total <- bind_rows(missing_table_data, any_row)
print(
  missing_table_with_total %>%
    select(characteristic, n_missing, pct_missing),
  n = Inf
)
missing_table_gt <- missing_table_with_total %>%
  select(
    Characteristic = characteristic,
    `Missing or unknown, n (%)` = `Missing/unknown, n (%)`
  ) %>%
  gt() %>%
  tab_header(
    title = md(
      "**Table S__. Missing or unknown covariate data among otherwise eligible infants.**"
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Characteristic,
      rows = Characteristic == "Any required covariate"
    )
  ) %>%
  tab_source_note(
    source_note = md(
      paste0(
        "Missingness was assessed among 47,868 infants meeting eligibility ",
        "criteria before complete-case exclusion. Missing or unknown values ",
        "were defined using the same covariate-completeness rules applied in ",
        "the Cox-model analysis. Individual infants could have missing or ",
        "unknown data for more than one covariate; therefore, ",
        "variable-specific counts need not sum to 218."
      )
    )
  )

missing_table_gt

gtsave(
  missing_table_gt,
  "table_missing_covariates_exact.docx"
)

write.csv(
  missing_table_with_total %>%
    select(characteristic, n_missing, pct_missing),
  "table_missing_covariates_exact.csv",
  row.names = FALSE
)

write.csv(
  discordant_infants,
  "missingness_validation_discordant_infants.csv",
  row.names = FALSE
)
