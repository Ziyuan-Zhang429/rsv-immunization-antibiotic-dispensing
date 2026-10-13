

stopifnot(
  exists("cohort_primary"), exists("model_data_primary"),
  exists("model_data_fullyear"), exists("model_data_out_of_season")
)

EXPOSURE_GROUP_LEVELS <- c(
  "Neither", "Maternal RSV vaccine only", "Nirsevimab only", "Both"
)

EXPOSURE_GROUP_DISPLAY <- c(
  "Neither"                   = "Neither",
  "Maternal RSV vaccine only" = "Maternal RSV vaccine",
  "Nirsevimab only"           = "Nirsevimab",
  "Both"                      = "Both"
)

canonical_exposure_group <- function(df) {
  stopifnot(all(c("rsv_valid", "mab_history") %in% names(df)))
  df %>%
    mutate(
      exposure_group = case_when(
        rsv_valid &  mab_history ~ "Both",
        rsv_valid                ~ "Maternal RSV vaccine only",
        mab_history              ~ "Nirsevimab only",
        TRUE                     ~ "Neither"
      ),
      exposure_group = factor(exposure_group, levels = EXPOSURE_GROUP_LEVELS)
    )
}

model_infant_ids <- function(md_list, outcome = "all_cause") {
  md_list[[outcome]] %>% distinct(infant_study_id)
}

complete_case_flag <- clean_and_restrict(cp_primary$all_cause, restrict_to_season = FALSE) %>%
  distinct(infant_study_id) %>%
  mutate(has_complete_covariates = TRUE)

analytic_infants <- cohort_primary %>%
  canonical_exposure_group() %>%
  left_join(complete_case_flag, by = "infant_study_id") %>%
  mutate(
    has_complete_covariates = replace_na(has_complete_covariates, FALSE),
    in_primary_model  = infant_study_id %in% model_infant_ids(model_data_primary)$infant_study_id,
    in_fullyear_model = infant_study_id %in% model_infant_ids(model_data_fullyear)$infant_study_id,
    in_oos_model      = infant_study_id %in% model_infant_ids(model_data_out_of_season)$infant_study_id
  )

stopifnot(

  all(analytic_infants$in_fullyear_model <= analytic_infants$has_complete_covariates)
)

stopifnot(
  n_distinct(analytic_infants$infant_study_id) == nrow(analytic_infants),
  !any(is.na(analytic_infants$exposure_group)),
  nrow(analytic_infants) == nrow(cohort_primary),

  sum(analytic_infants$exposure_group == "Both") == 0
)

analytic_infants <- analytic_infants %>%
  mutate(exposure_group = droplevels(exposure_group))

canonical_counts <- function(subset_flag = NULL, data = analytic_infants) {
  d <- if (is.null(subset_flag)) data else data %>% filter(.data[[subset_flag]])
  d %>%
    count(exposure_group, name = "n_infants_cohort") %>%
    complete(exposure_group, fill = list(n_infants_cohort = 0L)) %>%
    mutate(
      total = sum(n_infants_cohort),
      pct = 100 * n_infants_cohort / total,
      exposure = unname(EXPOSURE_GROUP_DISPLAY[as.character(exposure_group)])
    )
}

CANONICAL_N <- list(
  cohort         = canonical_counts(),
  complete_case  = canonical_counts("has_complete_covariates"),
  primary_model  = canonical_counts("in_primary_model"),
  fullyear_model = canonical_counts("in_fullyear_model"),
  oos_model      = canonical_counts("in_oos_model")
)
print(imap_dfr(CANONICAL_N, ~ .x %>% mutate(population = .y, .before = 1)) %>%
        select(population, exposure, n_infants_cohort, total, pct), n = Inf)

DAYS_PER_YEAR <- 365.25

format_exposure_arm <- function(x) {
  x <- sub("^exposure_arm", "", x)
  dplyr::case_when(
    x == "d_neither"           ~ "Neither",
    x == "a_maternal_vax_only" ~ "Maternal RSV vaccine",
    x == "b_nirsevimab_only"   ~ "Nirsevimab",
    x == "c_both"              ~ "Both",
    TRUE ~ NA_character_
  )
}

arm_counts <- function(md, label = NA_character_) {
  stopifnot(all(c("infant_study_id", "exposure_arm", "tstart", "tstop", "event") %in% names(md)))
  md %>%
    group_by(exposure_arm) %>%
    summarise(
      n_infants_contributing = n_distinct(infant_study_id),
      person_days = sum(tstop - tstart),
      n_events = sum(event),
      .groups = "drop"
    ) %>%
    mutate(
      person_years = person_days / DAYS_PER_YEAR,
      events_per_100_py = 100 * n_events / person_years,
      exposure = paste0("exposure_arm", exposure_arm),
      exposure_label = format_exposure_arm(as.character(exposure_arm)),
      model_label = label
    )
}

exclusive_counts <- function(md) {
  md %>%
    distinct(infant_study_id) %>%
    inner_join(analytic_infants %>% select(infant_study_id, exposure_group),
               by = "infant_study_id") %>%
    count(exposure_group, name = "n_infants_exclusive") %>%
    complete(exposure_group, fill = list(n_infants_exclusive = 0L)) %>%
    mutate(exposure_label = unname(EXPOSURE_GROUP_DISPLAY[as.character(exposure_group)]))
}

check_arm_counts <- function(md, label) {
  ac <- arm_counts(md, label)
  ex <- exclusive_counts(md)
  n_total <- n_distinct(md$infant_study_id)
  straddlers <- md %>%
    group_by(infant_study_id) %>%
    summarise(n_arms = n_distinct(exposure_arm), .groups = "drop") %>%
    filter(n_arms > 1) %>%
    nrow()
  excess <- sum(ac$n_infants_contributing) - n_total
  stopifnot(
    excess == straddlers,
    sum(ex$n_infants_exclusive) == n_total
  )
  invisible(list(contributing = ac, exclusive = ex, n_total = n_total,
                 excess = excess, straddlers = straddlers))
}

FIG_DENOM <- "infants"

format_events_denom <- function(n_events, n_infants_contributing, person_years,
                                denom = FIG_DENOM) {
  switch(
    denom,
    person_years = sprintf("%s / %s", scales::comma(n_events),
                           scales::comma(round(person_years))),
    infants      = sprintf("%s / %s", scales::comma(n_events),
                           scales::comma(n_infants_contributing)),
    stop("Unknown FIG_DENOM: ", denom)
  )
}

fig_denom_header <- function(denom = FIG_DENOM) {
  switch(denom,
         person_years = "Events / person-yr",
         infants      = "Events / infants",
         stop("Unknown FIG_DENOM: ", denom))
}

ARM_COUNT_NOTE <- paste0(
  "Counts are model-specific. Nirsevimab exposure is time-varying, so an ",
  "infant contributes pre-dose person-time to the neither category and ",
  "post-dose person-time to the nirsevimab category; infant counts by ",
  "exposure category are therefore not mutually exclusive and do not sum to ",
  "the cohort total. Events are attributed to the exposure state in effect ",
  "on the event day. Mutually exclusive infant counts are given in Table 1 ",
  "and Table S2."
)

fig_denom_note <- function(denom = FIG_DENOM) {
  if (denom == "person_years") {
    paste0("Person-years are accrued within each exposure state. ", ARM_COUNT_NOTE)
  } else {
    ARM_COUNT_NOTE
  }
}

mab_date_missing_n <- cohort_primary %>%
  filter(mab_yn == 1, is.na(mab_days)) %>%
  nrow()
invisible(check_arm_counts(model_data_primary$all_cause,       "primary (RSV season)"))
invisible(check_arm_counts(model_data_fullyear$all_cause,      "full year"))
invisible(check_arm_counts(model_data_out_of_season$all_cause, "out of season"))
