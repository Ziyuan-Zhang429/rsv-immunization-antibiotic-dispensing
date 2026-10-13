

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(readr)
  library(purrr); library(broom); library(survival); library(lubridate)
})

CANONICAL_SCRIPT <- 'R/primary/cox_primary.R'
HOSP_SAS <- 'kpsc_infant_hosp_10062026.sas7bdat'
LINK_FILE <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
OUT_DIR <- 'new_hospital_admissions_outputs'

MERGE_NEXT_DAY_ADJACENT <- FALSE

BIRTH_STAY_START_DAYS <- c(0, 1)

for (f in c(CANONICAL_SCRIPT, HOSP_SAS, LINK_FILE)) {
  if (!file.exists(f)) stop('Required file not found: ', f)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, name) { write_csv(x, file.path(OUT_DIR, name), na = ''); invisible(x) }
chr_id <- function(x) trimws(as.character(x))

source(CANONICAL_SCRIPT)
needed <- c('cohort_primary','blocks_primary','clean_and_restrict',
            'fit_ag','ag_formula','FOLLOWUP_MAX','EXPOSURE_LEVELS',
            'COVARIATE_TERMS')
for (nm in needed) if (!exists(nm, inherits = TRUE)) stop('Missing canonical object: ', nm)
stopifnot(FOLLOWUP_MAX == 180)

link <- read_csv(LINK_FILE, show_col_types = FALSE,
                 col_types = cols(.default = col_character())) %>%
  mutate(across(c(infant_study_id_0113, infant_study_id_0616), chr_id))
valid_link <- link %>%
  filter(has_new_id == 'TRUE', has_hospital_source_row == 'TRUE',
         hospital_history_valid == 'TRUE') %>%
  select(infant_study_id_0113, infant_study_id_0616) %>%
  semi_join(cohort_primary %>% transmute(infant_study_id_0113 = chr_id(infant_study_id)),
            by = 'infant_study_id_0113')
if (!nrow(valid_link) || anyDuplicated(valid_link$infant_study_id_0113) ||
    anyDuplicated(valid_link$infant_study_id_0616)) stop('Invalid or empty official linkage')

h <- read_sas(HOSP_SAS) %>% transmute(
  infant_study_id_0616 = chr_id(infant_study_id),
  admit_days = as.numeric(admit_days), disch_days = as.numeric(disch_days)) %>%
  inner_join(valid_link, by = 'infant_study_id_0616')
invalid <- h %>% filter(is.na(admit_days) | is.na(disch_days) |
                        admit_days < 0 | disch_days < admit_days)
if (nrow(invalid)) stop('Invalid linked hospitalization intervals: ', nrow(invalid),
                        '. Resolve against official QC before fitting.')
h <- h %>% distinct(infant_study_id_0113, admit_days, disch_days)

merge_one <- function(z) {
  z <- arrange(z, admit_days, disch_days)
  a <- numeric(); b <- numeric()
  for (i in seq_len(nrow(z))) {
    gap <- if (MERGE_NEXT_DAY_ADJACENT) 1 else 0
    if (!length(b) || z$admit_days[i] > tail(b, 1) + gap) {
      a <- c(a, z$admit_days[i]); b <- c(b, z$disch_days[i])
    } else b[length(b)] <- max(b[length(b)], z$disch_days[i])
  }
  tibble(admit_days = a, disch_days = b, episode_number = seq_along(a))
}
episodes <- bind_rows(lapply(split(h, h$infant_study_id_0113), merge_one),
                      .id = 'infant_study_id_0113') %>%
  group_by(infant_study_id_0113) %>%
  arrange(admit_days, disch_days, .by_group = TRUE) %>%
  mutate(presumed_birth_stay = episode_number == 1L &
           first(admit_days) %in% BIRTH_STAY_START_DAYS,
         new_admission = !presumed_birth_stay & admit_days > 0 &
           admit_days <= FOLLOWUP_MAX) %>% ungroup()

events <- episodes %>% filter(new_admission) %>%
  transmute(infant_study_id = infant_study_id_0113,
            event_day = as.numeric(admit_days)) %>% distinct()

if (is.numeric(cohort_primary$infant_study_id)) {
  events$infant_study_id <- as.numeric(events$infant_study_id)
} else {
  events$infant_study_id <- as.character(events$infant_study_id)
}
cohort <- blocks_primary$cohort %>%
  filter(chr_id(infant_study_id) %in% valid_link$infant_study_id_0113)
if (!nrow(cohort)) stop('No canonical infants remain after linkage')

occupancy <- episodes %>%
  filter(admit_days < disch_days, admit_days < FOLLOWUP_MAX, disch_days > 0) %>%
  transmute(infant_study_id = infant_study_id_0113,
            start_day = pmax(0, admit_days),
            end_day = pmin(FOLLOWUP_MAX, disch_days)) %>%
  filter(end_day > start_day)
if (is.numeric(cohort$infant_study_id)) {
  occupancy$infant_study_id <- as.numeric(occupancy$infant_study_id)
} else occupancy$infant_study_id <- as.character(occupancy$infant_study_id)
occupancy_changes <- bind_rows(
  occupancy %>% transmute(infant_study_id, time = start_day, occupied = 1L),
  occupancy %>% transmute(infant_study_id, time = end_day, occupied = 0L)
) %>% arrange(infant_study_id, time)
if (anyDuplicated(occupancy_changes[c('infant_study_id', 'time')])) {
  stop('Conflicting occupancy boundaries: inspect same-day adjacent episodes.')
}

base_surv <- cohort %>% transmute(infant_study_id, tstart = 0, tstop = FOLLOWUP_MAX)
cp <- tmerge(cohort %>% select(infant_study_id), base_surv,
             id = infant_study_id, tstop = tstop)
cp <- tmerge(cp, events, id = infant_study_id, event = event(event_day))
cp <- tmerge(cp, blocks_primary$mab_states, id = infant_study_id,
             mab_state = tdc(time, mab_state))
utd_tmerge <- blocks_primary$utd %>%
  mutate(utd_status_chr = as.character(utd_status))
cp <- tmerge(cp, utd_tmerge, id = infant_study_id,
             utd_status = tdc(interval_start, utd_status_chr))
cp <- tmerge(cp, blocks_primary$season_cuts, id = infant_study_id,
             season_state = tdc(time, season_state))
if (nrow(occupancy_changes)) {
  cp <- tmerge(cp, occupancy_changes, id = infant_study_id,
               hospital_occupied = tdc(time, occupied))
} else cp$hospital_occupied <- 0L

cp <- cp %>%
  mutate(mab_state = replace_na(mab_state, 0L),
         mab_status = factor(mab_state, levels = 0:3,
                             labels = c('unexposed','buffer','exposed','ambiguous_waned')),
         utd_status = factor(if_else(is.na(utd_status), 'Up to date', utd_status),
                             levels = c('Up to date','Not up to date','Unknown')),
         season_state = replace_na(season_state, 0L),
         in_season = season_state == 1L,
         hospital_occupied = replace_na(hospital_occupied, 0L)) %>%
  left_join(cohort %>% transmute(
    infant_study_id, index_date,
    exposure_maternal = if_else(rsv_valid, 'maternal_vax', 'no_maternal_vax'),
    age_index, race_cat_infant, race_cat_mom, gestage_imp, birth_ct,
    apgar_1min, apgar_5min, GENDER,
    charlson, charlson_cat, charlson_Wt, charlsonwt_cat,
    bmi_prior, bmi_cat, smoking, median_income_cat, OUTPT, INPT, ED,
    flu_prior, tdap_prior, tot_covid, lb_prior, gbs_prior, tx_ever,
    serious_newborn_comorbidity, mom_study_id), by = 'infant_study_id') %>%
  mutate(exposure_arm = case_when(
    mab_status %in% c('buffer','ambiguous_waned') ~ NA_character_,
    exposure_maternal == 'maternal_vax' & mab_status == 'exposed' ~ 'c_both',
    exposure_maternal == 'maternal_vax' & mab_status == 'unexposed' ~ 'a_maternal_vax_only',
    exposure_maternal == 'no_maternal_vax' & mab_status == 'exposed' ~ 'b_nirsevimab_only',
    exposure_maternal == 'no_maternal_vax' & mab_status == 'unexposed' ~ 'd_neither',
    TRUE ~ NA_character_),
    exposure_arm = factor(exposure_arm, levels = EXPOSURE_LEVELS),
    calendar_date_start = index_date + tstart,
    calendar_month = floor_date(calendar_date_start, 'month'))

if (any(cp$event == 1L & cp$hospital_occupied == 1L)) {
  stop('Admission event in occupied interval; inspect occupancy boundaries')
}
raw_event_count <- nrow(events)
cp_event_count <- sum(cp$event)
if (cp_event_count != raw_event_count) {
  stop('Admission event mismatch: raw=', raw_event_count,
       '; tmerge=', cp_event_count)
}
cp_risk <- cp %>% filter(hospital_occupied == 0L)
if (sum(cp_risk$event) != raw_event_count) stop('Lost admission events after occupancy filter')
stopifnot(all(cp_risk$tstop > cp_risk$tstart))

periods <- c('full_year','RSV_season','non_RSV_season')
crude_formula <- Surv(tstart,tstop,event) ~ exposure_arm + strata(calendar_month)
fit_period <- function(period) {
  md <- clean_and_restrict(cp_risk,
    restrict_to_season = (period == 'RSV_season'),
    keep_out_of_season_only = (period == 'non_RSV_season'))
  if (!nrow(md)) stop('No at-risk intervals: ', period)
  counts <- md %>% group_by(exposure_arm) %>% summarise(
    infants_contributing_time = n_distinct(infant_study_id),
    infants_with_admission = n_distinct(infant_study_id[event == 1]),
    admissions = sum(event), person_days = sum(tstop - tstart),
    .groups = 'drop') %>%
    mutate(period = period, person_years = person_days / 365.25, .before = 1)
  model_fit <- function(formula, model) {
    fit <- tryCatch(fit_ag(md, formula = formula), error = identity)
    if (inherits(fit, 'error')) return(tibble(period = period, model = model,
      term = NA_character_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
      status = paste('ERROR:', conditionMessage(fit))))
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(grepl('^exposure_arm', term)) %>%
      transmute(period = period, model = model, term = term,
                HR = estimate, CI_low = conf.low, CI_high = conf.high,
                VE_pct = 100 * (1 - estimate),
                VE_CI_low = 100 * (1 - conf.high),
                VE_CI_high = 100 * (1 - conf.low), status = 'OK')
  }
  list(counts = counts, ve = bind_rows(
    model_fit(crude_formula, 'calendar_stratified_unadjusted'),
    model_fit(ag_formula, 'fully_adjusted')),
    summary = tibble(period = period, infants = n_distinct(md$infant_study_id),
                     admissions = sum(md$event),
                     person_days = sum(md$tstop - md$tstart)))
}
results <- setNames(lapply(periods, fit_period), periods)
counts <- bind_rows(map(results, 'counts'))
ve <- bind_rows(map(results, 've'))
period_summary <- bind_rows(map(results, 'summary'))

qc <- tibble(metric = c('canonical_primary_infants', 'valid_officially_linked_infants',
  'distinct_hospital_intervals', 'merged_hospital_episodes',
  'presumed_birth_stays', 'infants_first_admission_after_day1',
  'new_admissions_days1_to180', 'infants_with_new_admission',
  'events_in_tmerge', 'events_after_occupancy_filter',
  'at_risk_person_days_before_model_exclusions',
  'occupied_person_days_before_model_exclusions'),
  value = c(nrow(blocks_primary$cohort), nrow(cohort), nrow(h), nrow(episodes),
    sum(episodes$presumed_birth_stay),
    episodes %>% group_by(infant_study_id_0113) %>%
      summarise(first_day = min(admit_days), .groups = 'drop') %>%
      summarise(n = sum(first_day > 1)) %>% pull(n),
    nrow(events), n_distinct(events$infant_study_id),
    cp_event_count, sum(cp_risk$event),
    sum(cp_risk$tstop - cp_risk$tstart),
    sum(cp$tstop - cp$tstart) - sum(cp_risk$tstop - cp_risk$tstart)))
wcsv(qc, '01_hospitalization_QC.csv')
wcsv(episodes %>% count(presumed_birth_stay, new_admission, name = 'episodes'),
     '02_episode_type_counts.csv')
wcsv(period_summary, '03_period_model_counts.csv')
wcsv(counts, '04_admissions_person_time_by_arm.csv')
wcsv(ve, '05_VE_unadjusted_and_adjusted.csv')
term_map <- c(exposure_arma_maternal_vax_only = 'a_maternal_vax_only',
              exposure_armb_nirsevimab_only = 'b_nirsevimab_only')
main <- ve %>% filter(model == 'fully_adjusted', term %in% names(term_map)) %>%
  mutate(exposure_arm = unname(term_map[term])) %>%
  left_join(counts, by = c('period','exposure_arm'))
ref <- counts %>% filter(as.character(exposure_arm) == 'd_neither') %>%
  transmute(period, reference_infants = infants_contributing_time,
            reference_admissions = admissions, reference_person_days = person_days)
main <- main %>% left_join(ref, by = 'period') %>% transmute(
  period, exposure_arm, reference_infants, reference_admissions,
  reference_person_days, exposed_infants = infants_contributing_time,
  exposed_admissions = admissions, exposed_person_days = person_days,
  VE_pct, VE_CI_low, VE_CI_high, status)
wcsv(main, '06_manuscript_style_adjusted_estimates.csv')
invisible(NULL); print(qc, n = Inf)
invisible(NULL); print(main, n = Inf)
