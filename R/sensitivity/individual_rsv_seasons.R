

library(haven)
library(tidyverse)
library(survival)
library(broom)
library(ggplot2)
library(forcats)
library(stringr)

select <- dplyr::select
filter <- dplyr::filter

VERBOSE <- FALSE

infant_cov      <- read_sas("kpsc_infant_covariates_01132026.sas7bdat")
infant_outcomes <- read_sas("kpsc_infant_outcomes_01132026.sas7bdat")
mom_cov         <- read_sas("kpsc_mom_covariates_01132026.sas7bdat")

stopifnot(
  n_distinct(infant_cov$infant_study_id) == nrow(infant_cov),
  n_distinct(mom_cov$infant_study_id) == nrow(mom_cov),
  all(mom_cov$infant_study_id %in% infant_cov$infant_study_id),
  all(infant_cov$infant_study_id %in% mom_cov$infant_study_id),
  all(unique(infant_outcomes$infant_study_id) %in% infant_cov$infant_study_id)
)

FOLLOWUP_MAX        <- 180
GESTAGE_BIRTH_MIN   <- 224
GESTAGE_AT_VAX_MIN  <- 224
GESTAGE_AT_VAX_MAX  <- 252
MIN_DAYS_VAX_TO_BIRTH <- 14

MAB_MIN_DAYS_POST   <- 0
MAB_MAX_DAYS_POST   <- 181
SIBLING_MATCH_TOLERANCE_DAYS <- 60
UTD_INTERVAL_DAYS   <- 30

RSV_SEASONS <- tribble(
  ~season_id, ~start,                ~end,
  1,          as.Date("2023-10-01"), as.Date("2024-02-29"),
  2,          as.Date("2024-11-01"), as.Date("2025-03-31")
)

in_rsv_season <- function(dates) {
  purrr::map_lgl(dates, function(d) {
    if (is.na(d)) return(NA)
    any(d >= RSV_SEASONS$start & d <= RSV_SEASONS$end)
  })
}

EXPOSURE_LEVELS <- c("d_neither", "a_maternal_vax_only", "b_nirsevimab_only", "c_both")

WASHOUT_DAYS_GRID  <- c(0, 7, 14, 28)

MIN_EVENTS_TO_FIT <- 50

MIN_EVENTS_PER_ARM_BOTH <- 5

diagnosis_flags <- c(
  "asthma_allergy_yn", "bronchitis_yn", "fever_yn", "influenza_yn",
  "bacterial_infection_yn", "other_infection_yn", "sup_om_yn", "nsup_om",
  "pharyngitis_yn", "pneumonia_yn", "sinusitis_yn", "UTI_yn",
  "skin_cond_yn", "URI_yn"
)

diagnosis_interpretation <- tribble(
  ~diagnosis,               ~dx_group,
  "sup_om_yn",              "respiratory_plausible",
  "nsup_om",                "respiratory_plausible",
  "bronchitis_yn",          "respiratory_plausible",
  "pneumonia_yn",           "respiratory_plausible",
  "sinusitis_yn",           "respiratory_plausible",
  "URI_yn",                 "respiratory_plausible",
  "pharyngitis_yn",         "respiratory_plausible",
  "influenza_yn",           "respiratory_plausible",
  "UTI_yn",                 "nonrespiratory_control",
  "skin_cond_yn",           "nonrespiratory_control",
  "fever_yn",               "nonspecific",
  "bacterial_infection_yn", "nonspecific",
  "other_infection_yn",     "nonspecific",
  "asthma_allergy_yn",      "nonspecific"
)

COMPOSITE_DX_FLAGS <- list(
  lrti = c("bronchitis_yn", "pneumonia_yn"),
  urti = c("pharyngitis_yn", "sinusitis_yn", "URI_yn"),
  aom  = c("sup_om_yn"),
  ari_composite = c(
    "bronchitis_yn", "pneumonia_yn",
    "pharyngitis_yn", "sinusitis_yn", "URI_yn",
    "sup_om_yn",
    "nsup_om",
    "influenza_yn"
  )
)

infant_cov_renamed <- infant_cov %>%
  rename(RACE_ETH_infant = RACE_ETH, race_cat_infant = race_cat)

mom_cov_renamed <- mom_cov %>%
  rename(RACE_ETH_mom = RACE_ETH, race_cat_mom = race_cat)

cohort_base <- infant_cov_renamed %>%
  left_join(mom_cov_renamed, by = c("infant_study_id", "mom_study_id", "index_date"))

cohort_step1 <- cohort_base %>% filter(age_index >= 16, age_index <= 49)
cohort_step2 <- cohort_step1 %>% filter(gestage_imp >= GESTAGE_BIRTH_MIN)

attrition <- tibble(
  step        = "0. Starting cohort",
  n_remaining = nrow(cohort_base),
  n_dropped   = NA_integer_
) %>%
  add_row(
    step        = "1. Maternal age 16-49",
    n_remaining = nrow(cohort_step1),
    n_dropped   = nrow(cohort_base) - nrow(cohort_step1)
  ) %>%
  add_row(
    step        = sprintf("2. Gestational age >= %d days (32 weeks)", GESTAGE_BIRTH_MIN),
    n_remaining = nrow(cohort_step2),
    n_dropped   = nrow(cohort_step1) - nrow(cohort_step2)
  ) %>%
  add_row(step = "3a. Continuous membership [ASSUMED]",
          n_remaining = nrow(cohort_step2), n_dropped = 0L) %>%
  add_row(step = "3b. Infant enrollment within 60 days of birth [ASSUMED]",
          n_remaining = nrow(cohort_step2), n_dropped = 0L)

eligible_cohort <- cohort_step2
print(attrition)
multiple_births <- mom_cov %>%
  group_by(mom_study_id, index_date) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  left_join(infant_cov %>% select(infant_study_id, birth_ct), by = "infant_study_id")
infant_outcomes_dedup <- infant_outcomes %>%
  distinct(infant_study_id, abx_days, abx_generic, abx_route, .keep_all = TRUE)

analytic_events <- eligible_cohort %>%
  select(infant_study_id, index_date) %>%
  left_join(infant_outcomes_dedup, by = c("infant_study_id", "index_date")) %>%
  mutate(
    exclude_route_primary = case_when(
      is.na(abx_generic) ~ NA,
      abx_route %in% c("Ophthalmic", "Otic", "Topical") ~ TRUE,
      abx_route == "" ~ NA,
      TRUE ~ FALSE
    ),

    exclude_route_sensitivity = case_when(
      is.na(abx_generic) ~ NA,
      abx_route == "Ophthalmic" & abx_days <= 1 ~ TRUE,
      abx_route %in% c("Otic", "Topical") ~ TRUE,
      abx_route == "" ~ NA,
      TRUE ~ FALSE
    ),
    analyzable_event_primary = case_when(
      is.na(abx_generic) ~ NA,
      is.na(exclude_route_primary) ~ NA,
      exclude_route_primary ~ FALSE,
      TRUE ~ TRUE
    ),
    analyzable_event_sensitivity = case_when(
      is.na(abx_generic) ~ NA,
      is.na(exclude_route_sensitivity) ~ NA,
      exclude_route_sensitivity ~ FALSE,
      TRUE ~ TRUE
    )
  )
print(analytic_events %>% count(exclude_route_primary))
print(analytic_events %>% count(exclude_route_sensitivity))

build_outcome_events <- function(analyzable_col) {
  base <- analytic_events %>%
    filter(.data[[analyzable_col]] == TRUE, abx_days > 0)

  composites <- imap(COMPOSITE_DX_FLAGS, function(flags, outcome_name) {
    base %>%
      filter(if_any(all_of(flags), ~ .x == 1)) %>%
      distinct(infant_study_id, abx_days)
  })

  c(
    list(all_cause = base %>% distinct(infant_study_id, abx_days)),
    composites
  )
}

events_primary     <- build_outcome_events("analyzable_event_primary")
events_sensitivity <- build_outcome_events("analyzable_event_sensitivity")
walk2(events_primary, names(events_primary), ~ invisible(NULL))
walk2(events_sensitivity, names(events_sensitivity), ~ invisible(NULL))
composite_overlap_check <- bind_rows(
  events_primary$lrti          %>% mutate(dx_group = "lrti"),
  events_primary$urti          %>% mutate(dx_group = "urti"),
  events_primary$aom           %>% mutate(dx_group = "aom")
) %>%
  count(infant_study_id, abx_days, name = "n_composites") %>%
  count(n_composites)
print(composite_overlap_check)

build_diagnosis_events <- function(flag_col, analyzable_col = "analyzable_event_primary") {
  analytic_events %>%
    filter(.data[[analyzable_col]] == TRUE, abx_days > 0, .data[[flag_col]] == 1) %>%
    distinct(infant_study_id, abx_days)
}

build_class_events <- function(class_name, analyzable_col = "analyzable_event_primary") {
  analytic_events %>%
    filter(.data[[analyzable_col]] == TRUE, abx_days > 0, abx_class == class_name) %>%
    distinct(infant_study_id, abx_days)
}

n_all_cause_events <- analytic_events %>%
  filter(analyzable_event_primary == TRUE, abx_days > 0) %>%
  nrow()

diagnosis_mix_primary <- analytic_events %>%
  filter(analyzable_event_primary == TRUE, abx_days > 0) %>%
  summarise(across(all_of(diagnosis_flags), ~ sum(.x == 1, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "diagnosis", values_to = "n_events") %>%
  mutate(pct_of_all_cause_events = 100 * n_events / n_all_cause_events) %>%
  left_join(diagnosis_interpretation, by = "diagnosis") %>%
  arrange(desc(n_events))
print(diagnosis_mix_primary, n = Inf)

no_flag_events <- analytic_events %>%
  filter(analyzable_event_primary == TRUE, abx_days > 0) %>%
  filter(if_all(all_of(diagnosis_flags), ~ is.na(.x) | .x == 0))
print(no_flag_events %>% count(abx_class, sort = TRUE), n = Inf)
print(no_flag_events %>% count(abx_route, sort = TRUE), n = Inf)

apply_washout <- function(events_df, washout_days) {
  events_df %>%
    arrange(infant_study_id, abx_days) %>%
    group_by(infant_study_id) %>%
    mutate(new_episode = c(TRUE, diff(abx_days) > washout_days)) %>%
    ungroup() %>%
    filter(new_episode) %>%
    select(-new_episode)
}

eligible_cohort <- eligible_cohort %>%
  mutate(
    implied_vax_date = index_date - RSV_days,
    gestage_at_vax   = gestage_imp - RSV_days,
    days_vax_to_birth = gestage_imp - gestage_at_vax,
    rsv_history      = rsv_prior == 1 & !is.na(RSV_days),
    rsv_valid = rsv_history &
      days_vax_to_birth >= MIN_DAYS_VAX_TO_BIRTH &
      gestage_at_vax >= GESTAGE_AT_VAX_MIN &
      gestage_at_vax <= GESTAGE_AT_VAX_MAX,
    rsv_ambiguous = rsv_history & !rsv_valid
  )
print(eligible_cohort %>%
        filter(rsv_history) %>%
        summarise(
          n_total                  = n(),
          n_gestage_window_only    = sum(gestage_at_vax >= GESTAGE_AT_VAX_MIN &
                                           gestage_at_vax <= GESTAGE_AT_VAX_MAX),
          n_min_days_only          = sum(days_vax_to_birth >= MIN_DAYS_VAX_TO_BIRTH),
          n_both_conditions        = sum(rsv_valid)
        ))
print(eligible_cohort %>% count(rsv_history, rsv_valid, rsv_ambiguous))

sibling_reclass_candidates <- eligible_cohort %>%
  filter(rsv_prior == 1, RSV_days < 0) %>%
  select(infant_study_id, mom_study_id, index_date, RSV_days, implied_vax_date) %>%
  inner_join(
    eligible_cohort %>% distinct(mom_study_id, sibling_index_date = index_date),
    by = "mom_study_id",
    relationship = "many-to-many"
  ) %>%
  filter(sibling_index_date != index_date) %>%
  mutate(days_from_sibling_index = abs(as.numeric(implied_vax_date - sibling_index_date))) %>%
  group_by(infant_study_id) %>%
  slice_min(days_from_sibling_index, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    sibling_match_confirmed =
      days_from_sibling_index <= SIBLING_MATCH_TOLERANCE_DAYS &
      sibling_index_date > index_date
  ) %>%
  filter(sibling_match_confirmed) %>%
  distinct(infant_study_id, sibling_match_confirmed)
eligible_cohort <- eligible_cohort %>%
  left_join(sibling_reclass_candidates, by = "infant_study_id") %>%
  mutate(
    sibling_match_confirmed = replace_na(sibling_match_confirmed, FALSE),
    rsv_valid     = if_else(sibling_match_confirmed, FALSE, rsv_valid),
    rsv_ambiguous = if_else(sibling_match_confirmed, FALSE, rsv_ambiguous)
  )
print(eligible_cohort %>%
        count(rsv_history, rsv_valid, rsv_ambiguous, sibling_match_confirmed))

residual_ambiguous_check <- eligible_cohort %>%
  filter(rsv_prior == 1, RSV_days < -365) %>%
  summarise(
    n              = n(),
    n_unresolved   = sum(rsv_ambiguous),
    pct_unresolved = 100 * mean(rsv_ambiguous)
  )
print(residual_ambiguous_check)

eligible_cohort <- eligible_cohort %>%
  mutate(mab_history = mab_yn == 1 & !is.na(mab_days))
print(eligible_cohort %>% count(mab_history))
build_mab_state_changes <- function(cohort, buffer_days, max_days = MAB_MAX_DAYS_POST) {
  windows <- cohort %>%
    filter(mab_history) %>%
    transmute(
      infant_study_id,
      dose_time    = pmax(pmin(mab_days, FOLLOWUP_MAX), 0),
      exposed_time = pmax(pmin(mab_days + buffer_days, FOLLOWUP_MAX), 0),
      waned_time   = pmax(pmin(mab_days + max_days, FOLLOWUP_MAX), 0)
    )

  changes <- if (buffer_days > 0) {
    bind_rows(
      windows %>% transmute(infant_study_id, time = dose_time,    mab_state = 1L),
      windows %>% transmute(infant_study_id, time = exposed_time, mab_state = 2L),
      windows %>% transmute(infant_study_id, time = waned_time,   mab_state = 3L)
    )
  } else {
    bind_rows(
      windows %>% transmute(infant_study_id, time = dose_time,  mab_state = 2L),
      windows %>% transmute(infant_study_id, time = waned_time, mab_state = 3L)
    )
  }

  changes %>%
    filter(time < FOLLOWUP_MAX) %>%
    arrange(infant_study_id, time)
}

eligible_cohort <- eligible_cohort %>%
  mutate(
    exposure_baseline_approx = case_when(
      rsv_ambiguous              ~ NA_character_,
      rsv_valid  &  mab_history  ~ "c_both",
      rsv_valid  & !mab_history  ~ "a_maternal_vax_only",
      !rsv_valid &  mab_history  ~ "b_nirsevimab_only",
      !rsv_valid & !mab_history  ~ "d_neither",
      TRUE                       ~ NA_character_
    )
  )
print(eligible_cohort %>%
        count(exposure_baseline_approx) %>%
        mutate(pct = 100 * n / sum(n)))

cohort_full <- eligible_cohort %>%
  filter(!rsv_ambiguous) %>%
  mutate(
    serious_newborn_comorbidity = as.integer(
      bronchopulm_dx == 1 | immunodef_dx == 1 | downsyn_dx == 1
    ),

    rsv_valid_anytime = rsv_history & !sibling_match_confirmed & RSV_days > 0,
    ever_both = rsv_valid & mab_history
  )
print(cohort_full %>%
        count(rsv_valid, rsv_valid_anytime) %>%
        mutate(pct = 100 * n / sum(n)))

ever_both_ids <- cohort_full %>% filter(ever_both) %>% pull(infant_study_id)

cohort_primary <- cohort_full %>% filter(!ever_both)
both_composition_check <- cohort_full %>%
  filter(rsv_valid) %>%
  group_by(dropped_ever_both = ever_both) %>%
  summarise(
    n                        = n(),
    median_days_vax_to_birth = median(days_vax_to_birth, na.rm = TRUE),
    pct_vax_within_21d_birth = 100 * mean(days_vax_to_birth < 21, na.rm = TRUE),
    pct_serious_comorbidity  = 100 * mean(serious_newborn_comorbidity, na.rm = TRUE),
    pct_preterm              = 100 * mean(gestage_imp < 259, na.rm = TRUE),
    mean_gestage             = mean(gestage_imp, na.rm = TRUE),
    mean_apgar_5min          = mean(apgar_5min, na.rm = TRUE),
    mean_OUTPT               = mean(OUTPT, na.rm = TRUE),
    mean_INPT                = mean(INPT, na.rm = TRUE),
    mean_ED                  = mean(ED, na.rm = TRUE),
    .groups = "drop"
  )
print(both_composition_check)
build_utd_at <- function(df, day_col) {
  df %>%
    mutate(
      utd_hepb = case_when(
        {{ day_col }} < 31 ~ 1L,
        {{ day_col }} < 91 ~ if_else(is.na(HBV_dose1_days), 0L,
                                     as.integer(HBV_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(HBV_dose2_30dmin_days), 0L,
                       as.integer(HBV_dose2_30dmin_days < {{ day_col }}))
      ),
      utd_rota = case_when(
        {{ day_col }} < 91 ~ 1L,
        {{ day_col }} < 151 ~ if_else(is.na(rota_dose1_days), 0L,
                                      as.integer(rota_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(rota_dose2_30dmin_days), 0L,
                       as.integer(rota_dose2_30dmin_days < {{ day_col }}))
      ),
      utd_dtap = case_when(
        {{ day_col }} < 91 ~ 1L,
        {{ day_col }} < 151 ~ if_else(is.na(DTaP_dose1_days), 0L,
                                      as.integer(DTaP_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(DTaP_dose2_30dmin_days), 0L,
                       as.integer(DTaP_dose2_30dmin_days < {{ day_col }}))
      ),
      utd_hib = case_when(
        {{ day_col }} < 91 ~ 1L,
        {{ day_col }} < 151 ~ if_else(is.na(Hib_dose1_days), 0L,
                                      as.integer(Hib_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(Hib_dose2_30dmin_days), 0L,
                       as.integer(Hib_dose2_30dmin_days < {{ day_col }}))
      ),
      utd_pcv = case_when(
        {{ day_col }} < 91 ~ 1L,
        {{ day_col }} < 151 ~ if_else(is.na(PCV_dose1_days), 0L,
                                      as.integer(PCV_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(PCV_dose2_30dmin_days), 0L,
                       as.integer(PCV_dose2_30dmin_days < {{ day_col }}))
      ),
      utd_ipv = case_when(
        {{ day_col }} < 91 ~ 1L,
        {{ day_col }} < 151 ~ if_else(is.na(IPV_dose1_days), 0L,
                                      as.integer(IPV_dose1_days < {{ day_col }})),
        TRUE ~ if_else(is.na(IPV_dose2_30dmin_days), 0L,
                       as.integer(IPV_dose2_30dmin_days < {{ day_col }}))
      ),

      utd_all_num = case_when(
        if_any(c(utd_hepb, utd_rota, utd_dtap, utd_hib, utd_pcv, utd_ipv),
               ~ .x == 0L) ~ 0L,
        if_any(c(utd_hepb, utd_rota, utd_dtap, utd_hib, utd_pcv, utd_ipv),
               ~ is.na(.x)) ~ NA_integer_,
        TRUE ~ 1L
      ),

      utd_status = factor(
        case_when(
          utd_all_num == 1L ~ "Up to date",
          utd_all_num == 0L ~ "Not up to date",
          TRUE              ~ "Unknown"
        ),
        levels = c("Up to date", "Not up to date", "Unknown")
      )
    ) %>%
    select(-utd_all_num)
}

build_season_cutpoints <- function(cohort) {
  cohort %>%
    select(infant_study_id, index_date) %>%
    crossing(RSV_SEASONS) %>%
    transmute(
      infant_study_id,

      on_time = as.numeric(start - index_date) - 1,

      off_time = as.numeric(end - index_date)
    ) %>%
    filter(
      on_time < FOLLOWUP_MAX,
      off_time > 0
    ) %>%
    mutate(
      on_time = pmax(on_time, 0),
      off_time = pmin(off_time, FOLLOWUP_MAX)
    ) %>%
    {
      bind_rows(
        select(
          .,
          infant_study_id,
          time = on_time
        ) %>%
          mutate(season_state = 1L),

        select(
          .,
          infant_study_id,
          time = off_time
        ) %>%
          mutate(season_state = 0L)
      )
    } %>%
    arrange(
      infant_study_id,
      time
    )
}

make_analysis_blocks <- function(cohort, buffer_days = MAB_MIN_DAYS_POST,
                                 label = "unnamed") {
  list(
    label       = label,
    cohort      = cohort,
    base_surv   = cohort %>% transmute(infant_study_id, tstart = 0, tstop = FOLLOWUP_MAX),
    mab_states  = build_mab_state_changes(cohort, buffer_days = buffer_days),
    season_cuts = build_season_cutpoints(cohort),
    utd = cohort %>%
      select(infant_study_id, starts_with("HBV_"), starts_with("rota_"),
             starts_with("DTaP_"), starts_with("Hib_"), starts_with("PCV_"),
             starts_with("IPV_")) %>%
      crossing(interval_start = seq(0, FOLLOWUP_MAX, by = UTD_INTERVAL_DAYS)) %>%
      build_utd_at(interval_start) %>%
      select(infant_study_id, interval_start, utd_status)
  )
}
print(
  cohort_primary %>%
    select(infant_study_id, starts_with("HBV_"), starts_with("rota_"),
           starts_with("DTaP_"), starts_with("Hib_"), starts_with("PCV_"),
           starts_with("IPV_")) %>%
    mutate(interval_start = 150) %>%
    build_utd_at(interval_start) %>%
    count(utd_status) %>%
    mutate(pct = round(100 * n / sum(n), 2))
)

blocks_primary <- make_analysis_blocks(cohort_primary, label = "never_both")
blocks_full    <- make_analysis_blocks(cohort_full,    label = "both_retained")
build_cp <- function(outcome_events, blocks = blocks_primary,
                     maternal_exposure_flag = "rsv_valid") {

  cohort <- blocks$cohort

  events_tmerge <- outcome_events %>%
    semi_join(cohort, by = "infant_study_id") %>%

    filter(abx_days > 0, abx_days <= FOLLOWUP_MAX) %>%
    rename(event_day = abx_days)

  utd_tmerge <- blocks$utd %>%
    mutate(utd_status_chr = as.character(utd_status))

  cp <- tmerge(data1 = cohort %>% select(infant_study_id),
               data2 = blocks$base_surv, id = infant_study_id, tstop = tstop)
  cp <- tmerge(cp, events_tmerge,      id = infant_study_id, event = event(event_day))
  cp <- tmerge(cp, blocks$mab_states,  id = infant_study_id, mab_state = tdc(time, mab_state))
  cp <- tmerge(cp, utd_tmerge,         id = infant_study_id, utd_status = tdc(interval_start, utd_status_chr))
  cp <- tmerge(cp, blocks$season_cuts, id = infant_study_id, season_state = tdc(time, season_state))

  cp <- cp %>%
    mutate(
      mab_state  = replace_na(mab_state, 0L),
      mab_status = factor(mab_state, levels = 0:3,
                          labels = c("unexposed", "buffer", "exposed", "ambiguous_waned")),
      utd_status = factor(if_else(is.na(utd_status), "Up to date", utd_status),
                          levels = c("Up to date", "Not up to date", "Unknown")),
      season_state = replace_na(season_state, 0L),
      in_season    = season_state == 1L
    ) %>%
    left_join(
      cohort %>%
        transmute(
          infant_study_id, index_date,
          exposure_maternal = if_else(.data[[maternal_exposure_flag]],
                                      "maternal_vax", "no_maternal_vax"),
          age_index, race_cat_infant, race_cat_mom, gestage_imp, birth_ct,
          apgar_1min, apgar_5min, GENDER,
          charlson, charlson_cat, charlson_Wt, charlsonwt_cat,
          bmi_prior, bmi_cat, smoking,
          median_income_cat, OUTPT, INPT, ED, flu_prior, tdap_prior,
          tot_covid, lb_prior, gbs_prior, tx_ever,
          serious_newborn_comorbidity, mom_study_id
        ),
      by = "infant_study_id"
    ) %>%
    mutate(
      exposure_arm = case_when(
        mab_status %in% c("buffer", "ambiguous_waned") ~ NA_character_,
        exposure_maternal == "maternal_vax"    & mab_status == "exposed"   ~ "c_both",
        exposure_maternal == "maternal_vax"    & mab_status == "unexposed" ~ "a_maternal_vax_only",
        exposure_maternal == "no_maternal_vax" & mab_status == "exposed"   ~ "b_nirsevimab_only",
        exposure_maternal == "no_maternal_vax" & mab_status == "unexposed" ~ "d_neither",
        TRUE ~ NA_character_
      ),
      exposure_arm = factor(exposure_arm, levels = EXPOSURE_LEVELS),
      calendar_date_start = index_date + tstart,
      calendar_month = floor_date(calendar_date_start, "month")
    )

  stopifnot(all(cp$tstop > cp$tstart))
  cp
}

apply_unknown_exclusions <- function(df) {
  df %>%
    filter(
      GENDER != "",
      median_income_cat != "4: Unknown",
      smoking != "Unknown",
      bmi_cat != "5: Unknown"
    )
}

clean_and_restrict <- function(cp, restrict_to_season = TRUE,
                               keep_out_of_season_only = FALSE,
                               drop_both_intervals = FALSE) {
  out <- cp
  if (keep_out_of_season_only) {
    out <- out %>% filter(!in_season)
  } else if (restrict_to_season) {
    out <- out %>% filter(in_season)
  }
  if (drop_both_intervals) {
    out <- out %>% filter(exposure_arm != "c_both")
  }

  out <- out %>%
    filter(!is.na(exposure_arm)) %>%
    apply_unknown_exclusions() %>%
    mutate(
      exposure_arm      = droplevels(exposure_arm),
      GENDER            = droplevels(factor(GENDER)),
      median_income_cat = droplevels(factor(median_income_cat)),
      smoking           = droplevels(factor(smoking)),
      bmi_cat           = droplevels(factor(bmi_cat)),
      calendar_month    = droplevels(factor(calendar_month))
    )

  model_predictors <- all.vars(as.formula(paste("~", COVARIATE_TERMS)))
  n_before <- n_distinct(out$infant_study_id)
  out <- out[complete.cases(out[, model_predictors, drop = FALSE]), ]
  n_after <- n_distinct(out$infant_study_id)

  if (n_before != n_after) {
  }

  out
}

cp_primary <- imap(events_primary, ~ build_cp(.x))
print(cp_primary$all_cause %>% count(mab_status))
print(cp_primary$all_cause %>% count(exposure_arm))
stopifnot(sum(cp_primary$all_cause$exposure_arm == "c_both", na.rm = TRUE) == 0)

COVARIATE_TERMS <- paste(
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED + flu_prior + tdap_prior +",
  "tot_covid + gbs_prior + utd_status + lb_prior"
)

ag_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_month) +", COVARIATE_TERMS
))

fit_ag <- function(md, formula = ag_formula) {

  coxph(formula, data = md, cluster = mom_study_id, ties = "efron", model = TRUE)
}

tidy_ve <- function(model, label) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^exposure_arm", term), !grepl(":", term)) %>%
    transmute(
      outcome = label, exposure = term,
      HR = estimate, CI_low = conf.low, CI_high = conf.high,
      VE_pct     = 100 * (1 - estimate),
      VE_CI_low  = 100 * (1 - conf.high),
      VE_CI_high = 100 * (1 - conf.low),
      p.value
    )
}

empty_ve_row <- function(outcome_label, status, n_events_raw = NA_integer_,
                         n_events_model = NA_integer_) {
  tibble(
    outcome = outcome_label, exposure = NA_character_, HR = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, VE_pct = NA_real_,
    VE_CI_low = NA_real_, VE_CI_high = NA_real_, p.value = NA_real_,
    n_events_raw = n_events_raw, n_events_model = n_events_model,
    fit_status = status
  )
}

report_ve_counts <- function(md_list, label_prefix = "") {
  imap_dfr(md_list, ~ tibble(
    outcome = .y,
    rows    = nrow(.x),
    infants = n_distinct(.x$infant_study_id),
    events  = sum(.x$event)
  )) %>%
    mutate(set = label_prefix, .before = 1)
}

R3_RSV_WINDOWS <- tribble(
  ~rsv_season, ~window_start,          ~window_end,
  "2023-24",   as.Date("2023-10-01"), as.Date("2024-02-29"),
  "2024-25",   as.Date("2024-11-01"), as.Date("2025-03-31")
)

restrict_calendar_window <- function(cp, start_date, end_date) {
  out <- cp %>%
    mutate(
      .win_lo = as.numeric(start_date - index_date) - 1,
      .win_hi = as.numeric(end_date   - index_date),
      .new_tstart = pmax(tstart, .win_lo),
      .new_tstop  = pmin(tstop,  .win_hi),
      .new_event = as.integer(event == 1L & tstop <= .win_hi & tstop > .win_lo)
    ) %>%
    filter(.new_tstop > .new_tstart) %>%
    mutate(
      tstart = .new_tstart,
      tstop  = .new_tstop,
      event  = .new_event,
      calendar_date_start = index_date + tstart + 1,
      calendar_month = floor_date(calendar_date_start, "month")
    ) %>%
    select(-.win_lo, -.win_hi, -.new_tstart, -.new_tstop, -.new_event)

  clean_and_restrict(out, restrict_to_season = FALSE)
}

run_r3_rsv_season <- function(rsv_season, window_start, window_end) {
  md <- restrict_calendar_window(cp_primary$all_cause, window_start, window_end)
  fit <- fit_ag(md, ag_formula)

  counts <- md %>%
    group_by(exposure_arm) %>%
    summarise(
      infants_contributing = n_distinct(infant_study_id),
      person_days = sum(tstop - tstart),
      antibiotic_events = sum(event),
      .groups = "drop"
    ) %>%
    mutate(
      exposure = case_when(
        exposure_arm == "d_neither" ~ "Neither",
        exposure_arm == "a_maternal_vax_only" ~ "Maternal RSVpreF",
        exposure_arm == "b_nirsevimab_only" ~ "Nirsevimab",
        TRUE ~ as.character(exposure_arm)
      )
    )

  est <- tidy_ve(fit, "all_cause") %>%
    mutate(
      exposure = case_when(
        grepl("a_maternal_vax_only", exposure) ~ "Maternal RSVpreF",
        grepl("b_nirsevimab_only", exposure) ~ "Nirsevimab",
        TRUE ~ exposure
      ),
      VE_95CI = sprintf("%.1f (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high)
    ) %>%
    select(exposure, VE_95CI)

  est <- bind_rows(
    tibble(exposure = "Neither", VE_95CI = "Reference"),
    est
  )

  counts %>%
    select(-exposure_arm) %>%
    left_join(est, by = "exposure") %>%
    mutate(
      rsv_season = rsv_season,
      window_start = window_start,
      window_end = window_end,
      .before = 1
    ) %>%
    arrange(factor(exposure,
                   levels = c("Neither", "Maternal RSVpreF", "Nirsevimab")))
}

r3_rsv_table <- pmap_dfr(
  R3_RSV_WINDOWS,
  ~ run_r3_rsv_season(..1, ..2, ..3)
) %>%
  mutate(
    rsv_season = factor(rsv_season, levels = c("2023-24", "2024-25")),
    exposure = factor(exposure,
                      levels = c("Neither", "Maternal RSVpreF", "Nirsevimab"))
  ) %>%
  arrange(rsv_season, exposure) %>%
  mutate(across(c(rsv_season, exposure), as.character))

r3_rsv_manuscript_table <- r3_rsv_table %>%
  transmute(
    `RSV season` = rsv_season,
    Exposure = exposure,
    `Infants contributing person-time` = infants_contributing,
    `Antibiotic events` = antibiotic_events,
    `VE, % (95% CI)` = VE_95CI
  )

figure3_rsv_expected <- tribble(
  ~exposure,          ~expected_events, ~expected_infants,
  "Neither",                     862L,            21685L,
  "Maternal RSVpreF",            293L,             6695L,
  "Nirsevimab",                  632L,            17677L
)

r3_reconciliation <- r3_rsv_table %>%
  group_by(exposure) %>%
  summarise(
    season_specific_events_sum = sum(antibiotic_events),
    season_specific_infants_sum = sum(infants_contributing),
    .groups = "drop"
  ) %>%
  left_join(figure3_rsv_expected, by = "exposure") %>%
  mutate(
    events_match = season_specific_events_sum == expected_events,
    infants_match = season_specific_infants_sum == expected_infants
  )
print(r3_rsv_manuscript_table, n = Inf)
print(r3_reconciliation, n = Inf)

stopifnot(
  all(r3_reconciliation$events_match),
  all(r3_reconciliation$infants_match)
)

write.csv(
  r3_rsv_manuscript_table,
  "reviewer3_RSV_seasons_only_TABLE.csv",
  row.names = FALSE
)

write.csv(
  r3_reconciliation,
  "reviewer3_RSV_seasons_reconciliation_CHECK.csv",
  row.names = FALSE
)
