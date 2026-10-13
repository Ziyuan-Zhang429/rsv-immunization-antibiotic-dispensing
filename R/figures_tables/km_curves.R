

library(glmnet)
library(patchwork)
library(purrr)
library(tidyr)

WEIGHT_TRUNC_PCTILE     <- c(0.01, 0.99)
POSITIVITY_START        <- as.Date("2023-11-15")

TARGET_N_CONTROLS       <- 4
CALIPER_GRID_DAYS       <- c(0L, 3L, 7L)
MATCH_ANY_THRESHOLD     <- 0.95
MATCH_TARGET_THRESHOLD  <- 0.90
SAME_DAY_ABX_IS_PREINDEX <- TRUE

PS_FLOOR                <- 0.01
PS_FLOOR_FLAG_PCT       <- 1.0
NIR_WEIGHT_TRUNC_PCTILE <- c(0.01, 0.99)

FOLLOWUP_XMAX  <- FOLLOWUP_MAX

stopifnot(FOLLOWUP_MAX == 180, exists("UTD_INTERVAL_DAYS"), UTD_INTERVAL_DAYS == 30)
PLOT_FONT_SIZE <- 13
RISK_DAYS      <- seq(0, FOLLOWUP_XMAX, UTD_INTERVAL_DAYS)

figure_colors <- c(
  "Neither"              = "grey40",
  "Maternal RSV vaccine" = "steelblue4",
  "Nirsevimab"           = "darkorange3"
)

set.seed(20260809)

stopifnot(
  exists("apply_unknown_exclusions"), exists("analytic_infants"),
  exists("arm_counts"), exists("canonical_exposure_group")
)

truncate_to_first_event <- function(cp, id_col) {
  cp %>%
    arrange(.data[[id_col]], tstart) %>%
    group_by(.data[[id_col]]) %>%
    mutate(had_prior_event = cumsum(lag(event, default = 0L))) %>%
    filter(had_prior_event == 0L) %>%
    select(-had_prior_event) %>%
    ungroup()
}

km_primary_events <- events_primary$all_cause %>%
  filter(!is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX) %>%
  semi_join(cohort_primary %>% distinct(infant_study_id), by = "infant_study_id") %>%
  distinct(infant_study_id, abx_days)

stopifnot(
  !anyDuplicated(km_primary_events %>% select(infant_study_id, abx_days)),
  all(km_primary_events$abx_days > 0),
  all(km_primary_events$abx_days <= FOLLOWUP_MAX)
)

maternal_km_events <- km_primary_events

cp_all_cause_recurrent <- build_cp(maternal_km_events)
build_maternal_weights <- function(md, trunc_pctile = WEIGHT_TRUNC_PCTILE) {
  md <- md %>%
    filter(exposure_arm %in% c("d_neither", "a_maternal_vax_only")) %>%
    mutate(exposure_arm = droplevels(factor(exposure_arm)))

  cov_rhs <- as.formula(paste("~ ns(tstart, df=3) + ns(as.numeric(calendar_date_start), df=3) +", COVARIATE_TERMS))
  num_rhs <- ~ ns(tstart, df=3) + ns(as.numeric(calendar_date_start), df=3)

  predictors <- c("tstart", "calendar_date_start", all.vars(as.formula(paste("~", COVARIATE_TERMS))))
  n_na <- sum(!complete.cases(md[, predictors, drop = FALSE]))
  if (n_na > 0) {
    md <- md[complete.cases(md[, predictors, drop = FALSE]), ]
  }

  build_x <- function(f, d) model.matrix(f, data = d)[, -1, drop = FALSE]
  y <- md$exposure_arm

  fit_den <- cv.glmnet(build_x(cov_rhs, md), y, family = "binomial", alpha = 0, standardize = TRUE, nfolds = 5)
  fit_num <- cv.glmnet(build_x(num_rhs, md), y, family = "binomial", alpha = 0, standardize = TRUE, nfolds = 5)

  p_den_1 <- as.numeric(predict(fit_den, newx = build_x(cov_rhs, md), s = "lambda.1se", type = "response"))
  p_num_1 <- as.numeric(predict(fit_num, newx = build_x(num_rhs, md), s = "lambda.1se", type = "response"))

  second_level <- levels(y)[2]
  p_den <- ifelse(as.character(y) == second_level, p_den_1, 1 - p_den_1)
  p_num <- ifelse(as.character(y) == second_level, p_num_1, 1 - p_num_1)

  p_den <- pmin(pmax(p_den, PS_FLOOR), 1 - PS_FLOOR)
  p_num <- pmin(pmax(p_num, PS_FLOOR), 1 - PS_FLOOR)

  md <- md %>% mutate(sw = p_num / p_den)

  bounds <- quantile(md$sw, trunc_pctile, na.rm = TRUE)
  md <- md %>%
    mutate(
      sw_trunc = pmin(pmax(sw, bounds[1]), bounds[2]),
      exposure_label = recode(as.character(exposure_arm),
                              "d_neither" = "Neither",
                              "a_maternal_vax_only" = "Maternal RSV vaccine")
    )
  md
}

fit_maternal_km <- function(variant = c("in_season", "full_year", "out_of_season")) {
  variant <- match.arg(variant)

  md <- switch(variant,
               in_season     = clean_and_restrict(cp_all_cause_recurrent, restrict_to_season = TRUE),
               full_year     = clean_and_restrict(cp_all_cause_recurrent, restrict_to_season = FALSE),
               out_of_season = clean_and_restrict(cp_all_cause_recurrent, keep_out_of_season_only = TRUE)
  )

  md <- md %>%
    filter(index_date >= POSITIVITY_START) %>%
    filter(exposure_arm %in% c("d_neither", "a_maternal_vax_only"))

  md <- md %>% truncate_to_first_event(id_col = "infant_study_id")

  md_w <- build_maternal_weights(md)

  fit <- survfit(
    Surv(tstart, tstop, event) ~ exposure_label,
    data = md_w, weights = sw_trunc,
    id = infant_study_id, robust = TRUE
  )
  print(md_w %>%
          group_by(exposure_label) %>%
          summarise(rows = n(), infants = n_distinct(infant_study_id),
                    person_days = sum(tstop - tstart), events = sum(event), .groups = "drop"))

  list(data = md_w, fit = fit, variant = variant)
}

maternal_in_season     <- fit_maternal_km("in_season")
maternal_full_year     <- fit_maternal_km("full_year")
maternal_out_of_season <- fit_maternal_km("out_of_season")

stopifnot(
  maternal_full_year$data %>%
    filter(event == 1L) %>%
    count(infant_study_id) %>%
    summarise(ok = all(n == 1L)) %>%
    pull(ok)
)
old_style_in_season_events <- clean_and_restrict(
  build_cp(maternal_km_events %>% group_by(infant_study_id) %>%
             slice_min(abx_days, n = 1, with_ties = FALSE) %>% ungroup()),
  restrict_to_season = TRUE
) %>%
  filter(index_date >= POSITIVITY_START, exposure_arm %in% c("d_neither", "a_maternal_vax_only")) %>%
  summarise(events = sum(event)) %>% pull(events)
new_style_in_season_events <- sum(maternal_in_season$data$event)
control_candidates <- cohort_primary %>%
  filter(!rsv_valid) %>%
  transmute(
    control_id = infant_study_id,
    control_mom_id = mom_study_id,
    control_birth_date = index_date,
    control_mab_history = mab_history,
    control_mab_days = if_else(
      mab_history & !is.na(mab_days) & mab_days >= -2,
      pmax(mab_days, 0),
      mab_days
    )
  )

stopifnot(
  nrow(control_candidates) ==
    sum(analytic_infants$exposure_group %in% c("Neither", "Nirsevimab only")),
  sum(control_candidates$control_mab_history) ==
    sum(analytic_infants$exposure_group == "Nirsevimab only")
)

make_candidate_pairs <- function(recipients, controls, window_days) {
  recipients %>%
    tidyr::crossing(birth_offset = seq(-window_days, window_days)) %>%
    mutate(candidate_birth_date = treated_birth_date + birth_offset) %>%
    inner_join(controls, by = c("candidate_birth_date" = "control_birth_date"), relationship = "many-to-many") %>%
    filter(control_id != treated_id) %>%
    mutate(
      control_birth_date = candidate_birth_date,
      birth_diff_days = as.numeric(control_birth_date - treated_birth_date),
      control_index_age = as.numeric(km_index_date - control_birth_date),
      age_diff_days = control_index_age - index_age,
      caliper = window_days
    ) %>%
    filter(
      control_index_age >= 0,
      !control_mab_history | is.na(control_mab_days) | control_mab_days > control_index_age
    )
}

select_matched_controls <- function(recipients) {
  pairs <- map(CALIPER_GRID_DAYS, ~ make_candidate_pairs(recipients, control_candidates, .x))
  walk(pairs, ~ stopifnot(all(abs(.x$age_diff_days + .x$birth_diff_days) < 1e-8)))

  counts <- recipients %>%
    select(treated_id) %>%
    left_join(pairs[[1]] %>% count(treated_id, name = "n_c0"), by = "treated_id") %>%
    left_join(pairs[[2]] %>% count(treated_id, name = "n_c3"), by = "treated_id") %>%
    left_join(pairs[[3]] %>% count(treated_id, name = "n_c7"), by = "treated_id") %>%
    mutate(across(c(n_c0, n_c3, n_c7), ~ replace_na(.x, 0L))) %>%
    mutate(
      chosen_caliper = case_when(
        n_c0 >= TARGET_N_CONTROLS ~ 0L,
        n_c3 >= TARGET_N_CONTROLS ~ 3L,
        n_c7 >= TARGET_N_CONTROLS ~ 7L,
        n_c0 >= 1 ~ 0L,
        n_c3 >= 1 ~ 3L,
        n_c7 >= 1 ~ 7L,
        TRUE ~ NA_integer_
      ),
      n_available = case_when(
        chosen_caliper == 0L ~ n_c0,
        chosen_caliper == 3L ~ n_c3,
        chosen_caliper == 7L ~ n_c7,
        TRUE ~ 0L
      )
    )

  coverage <- counts %>%
    summarise(
      n_recipients = n(),
      pct_matched  = 100 * mean(!is.na(chosen_caliper)),
      pct_target   = 100 * mean(n_available >= TARGET_N_CONTROLS),
      pct_exact    = 100 * mean(chosen_caliper == 0L, na.rm = TRUE),
      pct_3day     = 100 * mean(chosen_caliper == 3L, na.rm = TRUE),
      pct_7day     = 100 * mean(chosen_caliper == 7L, na.rm = TRUE)
    )

  print(coverage)

  if (coverage$pct_matched < 100 * MATCH_ANY_THRESHOLD)
    warning("Fewer than the prespecified proportion of nirsevimab recipients could be matched.")
  if (coverage$pct_target < 100 * MATCH_TARGET_THRESHOLD)
    warning("Fewer than the prespecified proportion of recipients achieved the target control ratio.")

  candidates <- bind_rows(pairs) %>%
    inner_join(counts %>% select(treated_id, chosen_caliper), by = "treated_id") %>%
    filter(!is.na(chosen_caliper), caliper == chosen_caliper)

  set.seed(20260809)
  matched <- candidates %>%
    mutate(tie_break = runif(n())) %>%
    group_by(treated_id) %>%
    arrange(abs(age_diff_days), tie_break, .by_group = TRUE) %>%
    slice_head(n = TARGET_N_CONTROLS) %>%
    ungroup()

  list(matched = matched, counts = counts, coverage = coverage)
}

build_all_postindex_abx <- function(indexed_df, same_day_is_preindex = TRUE) {
  op <- if (same_day_is_preindex) `>` else `>=`

  indexed_df %>%
    select(episode_id, infant_study_id, index_age) %>%
    left_join(km_primary_events, by = "infant_study_id", relationship = "many-to-many") %>%
    filter(!is.na(abx_days), op(abx_days, index_age)) %>%
    transmute(episode_id, abx_day_since_index = abx_days - index_age) %>%
    distinct(episode_id, abx_day_since_index)
}

build_nirsevimab_indexed_unified <- function() {
  recipients <- cohort_primary %>%
    filter(
      !rsv_valid,
      mab_history,
      !is.na(mab_days),
      mab_days >= -2,
      mab_days < FOLLOWUP_MAX
    ) %>%
    mutate(
      mab_days_analysis = pmax(mab_days, 0)
    ) %>%
    transmute(
      treated_id = infant_study_id,
      treated_mom_id = mom_study_id,
      treated_birth_date = index_date,
      index_age = mab_days_analysis,
      km_index_date = index_date + mab_days_analysis
    )

  match_obj <- select_matched_controls(recipients)
  controls <- match_obj$matched

  exposed <- recipients %>%
    semi_join(controls %>% distinct(treated_id), by = "treated_id") %>%
    transmute(
      match_id = treated_id, infant_study_id = treated_id, mom_study_id = treated_mom_id,
      birth_date = treated_birth_date, km_index_date, index_age,
      exposure_km = "Nirsevimab", future_mab_day = NA_real_
    )

  unexposed <- controls %>%
    transmute(
      match_id = treated_id, infant_study_id = control_id, mom_study_id = control_mom_id,
      birth_date = control_birth_date, km_index_date, index_age = control_index_age,
      exposure_km = "Neither",
      future_mab_day = if_else(control_mab_history & control_mab_days > control_index_age,
                               control_mab_days, NA_real_)
    )

  dat <- bind_rows(exposed, unexposed) %>%
    mutate(exposure_km = factor(exposure_km, levels = c("Neither", "Nirsevimab")),
           episode_id = row_number())

  dat %>%
    mutate(
      days_to_age180    = FOLLOWUP_MAX - index_age,
      days_to_later_nir = if_else(exposure_km == "Neither" & !is.na(future_mab_day),
                                  future_mab_day - index_age, Inf),
      raw_max_followup  = pmin(days_to_age180, days_to_later_nir, na.rm = TRUE),
      censor_reason = if_else(
        exposure_km == "Neither" & is.finite(days_to_later_nir) &
          abs(days_to_later_nir - raw_max_followup) < 1e-8,
        "Control received nirsevimab", "Age 180 days"
      )
    ) %>%
    filter(raw_max_followup > 0)
}

add_nir_covariates <- function(dat) {
  covs <- dat %>%
    select(episode_id, infant_study_id, index_age) %>%
    left_join(
      cohort_primary %>%
        select(
          infant_study_id, GENDER, race_cat_infant, serious_newborn_comorbidity,
          apgar_1min, age_index, charlson_Wt, bmi_prior, bmi_cat, smoking,
          median_income_cat, OUTPT, INPT, ED, flu_prior, tdap_prior,
          tot_covid, gbs_prior, lb_prior,
          starts_with("HBV_"), starts_with("rota_"), starts_with("DTaP_"),
          starts_with("Hib_"), starts_with("PCV_"), starts_with("IPV_")
        ),
      by = "infant_study_id"
    ) %>%
    build_utd_at(index_age)

  prior_abx <- dat %>%
    select(episode_id, infant_study_id, index_age) %>%
    left_join(km_primary_events, by = "infant_study_id", relationship = "many-to-many") %>%
    group_by(episode_id) %>%
    summarise(
      n_prior_abx = sum(abx_days < index_age, na.rm = TRUE),
      any_prior_abx = as.integer(n_prior_abx > 0),
      .groups = "drop"
    )

  dat %>%
    left_join(
      covs %>%
        select(
          episode_id, GENDER, race_cat_infant, serious_newborn_comorbidity,
          apgar_1min, age_index, charlson_Wt, bmi_prior, bmi_cat, smoking,
          median_income_cat, OUTPT, INPT, ED, flu_prior, tdap_prior,
          tot_covid, gbs_prior, lb_prior, utd_status
        ),
      by = "episode_id"
    ) %>%
    left_join(prior_abx, by = "episode_id")
}

nir_ps_formula <- as.formula(paste(
  "I(exposure_km == 'Nirsevimab') ~",
  "GENDER + race_cat_infant + serious_newborn_comorbidity +",
  "apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +",
  "median_income_cat + OUTPT + INPT + ED + flu_prior + tdap_prior +",
  "tot_covid + gbs_prior + utd_status + lb_prior + any_prior_abx"
))

weight_nirsevimab_episode <- function(dat) {
  predictors <- all.vars(as.formula(paste("~", as.character(nir_ps_formula)[3])))
  complete <- complete.cases(dat[, predictors, drop = FALSE])

  if (sum(!complete) > 0)
  dat <- dat[complete, ]
  ps_fit <- glm(nir_ps_formula, data = dat, family = binomial())

  dat <- dat %>%
    mutate(
      ps_raw = predict(ps_fit, newdata = dat, type = "response"),
      ps_floored = ps_raw < PS_FLOOR | ps_raw > 1 - PS_FLOOR,
      ps = pmin(pmax(ps_raw, PS_FLOOR), 1 - PS_FLOOR)
    )

  pct_floored <- 100 * mean(dat$ps_floored)
  if (pct_floored > PS_FLOOR_FLAG_PCT)
    warning("Nirsevimab PS overlap failed prespecified flooring threshold.")

  p_treated <- mean(dat$exposure_km == "Nirsevimab")

  dat <- dat %>%
    mutate(
      sw_raw = if_else(exposure_km == "Nirsevimab", p_treated / ps, (1 - p_treated) / (1 - ps))
    )

  bounds <- quantile(dat$sw_raw, NIR_WEIGHT_TRUNC_PCTILE, na.rm = TRUE)

  dat <- dat %>%
    mutate(sw_trunc = pmin(pmax(sw_raw, bounds[1]), bounds[2]),
           exposure_label = as.character(exposure_km))
  dat
}

build_nir_recurrent_segments <- function(indexed_df) {

  abx <- build_all_postindex_abx(
    indexed_df,
    SAME_DAY_ABX_IS_PREINDEX
  ) %>%
    inner_join(
      indexed_df %>%
        select(
          episode_id,
          raw_max_followup
        ),
      by = "episode_id"
    ) %>%
    filter(
      abx_day_since_index > 0,
      abx_day_since_index <= raw_max_followup
    ) %>%
    distinct(
      episode_id,
      abx_day_since_index
    )

  abx_by_episode <- split(
    abx$abx_day_since_index,
    abx$episode_id
  )

  season_bounds_for_episode <- function(
    km_index_date,
    raw_max_followup
  ) {

    bounds <- c(
      as.numeric(
        RSV_SEASONS$start - km_index_date
      ) - 1,

      as.numeric(
        RSV_SEASONS$end - km_index_date
      )
    )

    sort(
      unique(
        bounds[
          bounds > 0 &
            bounds < raw_max_followup
        ]
      )
    )
  }

  indexed_df %>%
    mutate(
      .season_cuts = map2(
        km_index_date,
        raw_max_followup,
        season_bounds_for_episode
      ),

      .abx_cuts = map(
        episode_id,
        ~ {
          z <- abx_by_episode[[as.character(.x)]]

          if (is.null(z)) {
            numeric(0)
          } else {
            z
          }
        }
      ),

      .cuts = pmap(
        list(
          .season_cuts,
          .abx_cuts,
          raw_max_followup
        ),
        ~ sort(
          unique(
            c(
              0,
              ..1,
              ..2,
              ..3
            )
          )
        )
      )
    ) %>%
    select(
      -.season_cuts,
      -.abx_cuts
    ) %>%
    tidyr::unnest_longer(
      .cuts,
      values_to = "cutpoint"
    ) %>%
    group_by(episode_id) %>%
    arrange(
      episode_id,
      cutpoint
    ) %>%
    mutate(
      tstart = cutpoint,
      tstop = lead(cutpoint)
    ) %>%
    ungroup() %>%
    filter(
      !is.na(tstop),
      tstop > tstart
    ) %>%
    left_join(
      abx %>%
        mutate(event = 1L) %>%
        rename(
          tstop = abx_day_since_index
        ),
      by = c(
        "episode_id",
        "tstop"
      )
    ) %>%
    mutate(
      event = replace_na(event, 0L),

      seg_calendar_date = km_index_date + tstart + 1,

      seg_in_season = in_rsv_season(
        seg_calendar_date
      )
    ) %>%
    select(-cutpoint)
}

build_nirsevimab_base <- function() {
  indexed <- build_nirsevimab_indexed_unified() %>%
    add_nir_covariates() %>%
    apply_unknown_exclusions() %>%
    weight_nirsevimab_episode()

  build_nir_recurrent_segments(indexed)
}

nir_base <- build_nirsevimab_base()
print(
  nir_base %>%
    distinct(episode_id, exposure_label, censor_reason) %>%
    count(exposure_label, censor_reason) %>%
    group_by(exposure_label) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup(),
  n = Inf
)
fit_nirsevimab_km <- function(variant = c("in_season", "full_year", "out_of_season"),
                              base = nir_base) {
  variant <- match.arg(variant)

  dat <- switch(variant,
                in_season     = base %>% filter(seg_in_season),
                out_of_season = base %>% filter(!seg_in_season),
                full_year     = base
  )

  dat <- dat %>% truncate_to_first_event(id_col = "episode_id")

  fit <- survfit(
    Surv(tstart, tstop, event) ~ exposure_label,
    data = dat, weights = sw_trunc,
    id = episode_id, robust = TRUE
  )
  print(dat %>%
          group_by(exposure_label) %>%
          summarise(segments = n(), episodes = n_distinct(episode_id),
                    person_days = sum(tstop - tstart), events = sum(event), .groups = "drop"))

  list(data = dat, fit = fit, variant = variant)
}

nir_in_season     <- fit_nirsevimab_km("in_season")
nir_full_year     <- fit_nirsevimab_km("full_year")
nir_out_of_season <- fit_nirsevimab_km("out_of_season")

extract_km_curve <- function(fit, dat, arms, xmax) {
  sf <- summary(fit, censored = FALSE)

  df <- tibble(
    exposure = str_remove(as.character(sf$strata), "^exposure_label="),
    time = sf$time, surv = sf$surv, lower = sf$lower, upper = sf$upper
  ) %>%
    filter(exposure %in% arms) %>%
    mutate(cuminc = 1 - surv, cuminc_lower = 1 - upper, cuminc_upper = 1 - lower)

  origin <- df %>%
    distinct(exposure) %>%
    mutate(time = 0, surv = 1, lower = 1, upper = 1,
           cuminc = 0, cuminc_lower = 0, cuminc_upper = 0)

  max_times <- dat %>%
    group_by(exposure = exposure_label) %>%
    summarise(max_time = max(tstop), .groups = "drop")

  end <- df %>%
    group_by(exposure) %>%
    slice_max(time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(max_times, by = "exposure") %>%
    mutate(time = pmax(time, max_time, xmax)) %>%
    select(-max_time)

  bind_rows(origin, df, end) %>% arrange(exposure, time)
}

get_cuminc_at_end <- function(fit_obj, arms) {

  curve <- extract_km_curve(
    fit_obj$fit,
    fit_obj$data,
    arms = arms,
    xmax = FOLLOWUP_XMAX
  )

  curve %>%
    filter(time <= FOLLOWUP_XMAX) %>%
    group_by(exposure) %>%
    slice_max(time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(exposure, cuminc)
}

get_risk_difference <- function(fit_obj, exposed_arm) {

  risks <- get_cuminc_at_end(
    fit_obj,
    arms = c("Neither", exposed_arm)
  )

  risk_neither <- risks %>%
    filter(exposure == "Neither") %>%
    pull(cuminc)

  risk_exposed <- risks %>%
    filter(exposure == exposed_arm) %>%
    pull(cuminc)

  tibble(
    exposure = exposed_arm,
    risk_neither_pct = 100 * risk_neither,
    risk_exposed_pct = 100 * risk_exposed,
    absolute_difference_pct = 100 * (risk_neither - risk_exposed)
  )
}

km_absolute_differences <- bind_rows(

  get_risk_difference(
    maternal_full_year,
    "Maternal RSV vaccine"
  ) %>%
    mutate(period = "Full year"),

  get_risk_difference(
    nir_full_year,
    "Nirsevimab"
  ) %>%
    mutate(period = "Full year"),

  get_risk_difference(
    maternal_in_season,
    "Maternal RSV vaccine"
  ) %>%
    mutate(period = "RSV season"),

  get_risk_difference(
    nir_in_season,
    "Nirsevimab"
  ) %>%
    mutate(period = "RSV season")

) %>%
  select(
    period,
    exposure,
    risk_neither_pct,
    risk_exposed_pct,
    absolute_difference_pct
  )

print(km_absolute_differences, n = Inf)

stopifnot(
  nir_full_year$data %>%
    filter(event == 1L) %>%
    count(episode_id) %>%
    summarise(ok = all(n == 1L)) %>%
    pull(ok)
)
old_style_first_postindex <- function(indexed_df, same_day_is_preindex = TRUE) {
  op <- if (same_day_is_preindex) `>` else `>=`
  indexed_df %>%
    select(episode_id, infant_study_id, index_age) %>%
    left_join(km_primary_events, by = "infant_study_id", relationship = "many-to-many") %>%
    filter(is.na(abx_days) | op(abx_days, index_age)) %>%
    group_by(episode_id) %>%
    summarise(first_abx_day = if (all(is.na(abx_days))) NA_real_ else min(abx_days, na.rm = TRUE), .groups = "drop")
}
old_indexed <- build_nirsevimab_indexed_unified() %>%
  left_join(old_style_first_postindex(., SAME_DAY_ABX_IS_PREINDEX), by = "episode_id") %>%
  mutate(
    days_to_abx = first_abx_day - index_age,
    event_old   = as.integer(!is.na(days_to_abx) & days_to_abx <= raw_max_followup),
    seg_calendar_date_of_event = km_index_date + days_to_abx,
    event_old_in_season = event_old == 1L & in_rsv_season(seg_calendar_date_of_event)
  )
old_style_in_season_events <- sum(old_indexed$event_old_in_season, na.rm = TRUE)
new_style_in_season_events <- sum(nir_in_season$data$event)
km_risk_table <- function(dat, days = RISK_DAYS) {

  map_dfr(days, function(d) {
    at_risk <- if (d == 0) {
      dat %>% filter(tstart <= 0, tstop > 0)
    } else {
      dat %>% filter(tstart < d, tstop >= d)
    }

    at_risk %>%
      distinct(across(any_of(c("infant_study_id", "episode_id"))), exposure_label) %>%
      count(exposure_label, name = "n_at_risk") %>%
      transmute(exposure = exposure_label, day = d, n_at_risk)
  })
}

plot_final_km <- function(df, panel_title, x_title, xmax, arms,
                          y_max = 0.12, y_breaks = seq(0, 0.12, 0.02)) {
  ggplot(df, aes(x = time, y = cuminc, color = exposure, fill = exposure)) +
    geom_ribbon(aes(ymin = cuminc_lower, ymax = cuminc_upper),
                alpha = 0.15, color = NA, outline.type = "upper", show.legend = FALSE) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = figure_colors, breaks = arms, name = NULL) +
    scale_fill_manual(values = figure_colors, breaks = arms, name = NULL) +
    scale_y_continuous(
      breaks = y_breaks, limits = c(0, y_max),
      labels = scales::label_number(accuracy = 1, scale = 100),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(breaks = seq(0, xmax, 30), expand = expansion(mult = c(0.01, 0.02))) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    labs(title = panel_title, x = x_title, y = "Cumulative incidence (%)") +
    theme_classic(base_size = PLOT_FONT_SIZE) +
    theme(
      text = element_text(size = PLOT_FONT_SIZE),
      plot.title = element_text(size = PLOT_FONT_SIZE, face = "bold"),
      axis.title = element_text(size = PLOT_FONT_SIZE),
      axis.text = element_text(size = PLOT_FONT_SIZE, color = "black"),
      legend.text = element_text(size = PLOT_FONT_SIZE),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      plot.title.position = "plot",
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()
    )
}

plot_km_with_risk <- function(km_plot, risk_table, arms, x_title, xmax, risk_days = RISK_DAYS) {
  risk_df <- risk_table %>%
    filter(day %in% risk_days, exposure %in% arms) %>%
    mutate(exposure = factor(exposure, levels = rev(arms)), label = scales::comma(n_at_risk))

  risk_plot <- ggplot(risk_df, aes(x = day, y = exposure)) +
    geom_text(aes(label = label), color = "black", size = PLOT_FONT_SIZE / .pt) +
    scale_x_continuous(breaks = risk_days, expand = expansion(mult = c(0.05, 0.02))) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(x = x_title, y = NULL, title = "Number at risk") +
    theme_classic(base_size = PLOT_FONT_SIZE) +
    theme(
      text = element_text(size = PLOT_FONT_SIZE),
      plot.title = element_text(size = PLOT_FONT_SIZE, face = "bold", hjust = 0),
      axis.title.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
      axis.text.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
      axis.text.y = element_text(size = PLOT_FONT_SIZE, color = "black", hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
      plot.margin = margin(t = 2, b = 5, l = 5, r = 15)
    ) +
    coord_cartesian(clip = "off")

  km_top <- km_plot +
    labs(x = NULL) +
    theme(
      axis.text.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
      axis.title.y = element_text(size = PLOT_FONT_SIZE, margin = margin(r = 4)),
      plot.margin = margin(t = 5, b = 0, l = 5, r = 15)
    ) +
    coord_cartesian(clip = "off")

  km_top / risk_plot + plot_layout(heights = c(4, 1))
}

auto_y_scale <- function(...) {
  curves <- list(...)
  y_max <- max(map_dbl(curves, ~ max(.x$cuminc_upper, na.rm = TRUE))) * 1.1
  y_max <- ceiling(y_max / 0.02) * 0.02
  list(y_max = y_max, y_breaks = seq(0, y_max - 0.02, 0.02))
}

build_period_stack <- function(full_year_obj, in_season_obj, x_title, arms,
                               top_label = "Full year",
                               bottom_label = "RSV season",
                               risk_days = RISK_DAYS) {
  top_curve    <- extract_km_curve(full_year_obj$fit, full_year_obj$data, arms = arms, xmax = FOLLOWUP_XMAX)
  bottom_curve <- extract_km_curve(in_season_obj$fit, in_season_obj$data, arms = arms, xmax = FOLLOWUP_XMAX)

  p_top    <- plot_final_km(top_curve,    top_label,    x_title, FOLLOWUP_XMAX, arms)
  p_bottom <- plot_final_km(bottom_curve, bottom_label, x_title, FOLLOWUP_XMAX, arms)

  style_curve_top <- function(p) {
    p + labs(x = NULL) +
      theme(
        axis.text.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
        axis.title.y = element_text(size = PLOT_FONT_SIZE, margin = margin(r = 4)),
        plot.margin = margin(t = 5, b = 0, l = 5, r = 15)
      ) +
      coord_cartesian(clip = "off")
  }

  build_risk_leaf <- function(fit_data, x_title) {
    risk_df <- km_risk_table(fit_data) %>%
      filter(day %in% risk_days, exposure %in% arms) %>%
      mutate(exposure = factor(exposure, levels = rev(arms)), label = scales::comma(n_at_risk))

    ggplot(risk_df, aes(x = day, y = exposure)) +
      geom_text(aes(label = label), color = "black", size = PLOT_FONT_SIZE / .pt) +
      scale_x_continuous(breaks = risk_days, expand = expansion(mult = c(0.05, 0.02))) +
      scale_y_discrete(expand = expansion(add = 0.6)) +
      labs(x = x_title, y = NULL, title = "Number at risk") +
      theme_classic(base_size = PLOT_FONT_SIZE) +
      theme(
        text = element_text(size = PLOT_FONT_SIZE),
        plot.title = element_text(size = PLOT_FONT_SIZE, face = "bold", hjust = 0),
        axis.title.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
        axis.text.x = element_text(size = PLOT_FONT_SIZE, color = "black"),
        axis.text.y = element_text(size = PLOT_FONT_SIZE, color = "black", hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(t = 2, b = 5, l = 5, r = 15)
      ) +
      coord_cartesian(clip = "off")
  }

  list(
    top_curve    = style_curve_top(p_top),
    top_risk     = build_risk_leaf(full_year_obj$data, x_title),
    bottom_curve = style_curve_top(p_bottom),
    bottom_risk  = build_risk_leaf(in_season_obj$data, x_title)
  )
}

mat_leaves <- build_period_stack(
  maternal_full_year, maternal_in_season,
  x_title = "Days since birth",
  arms = c("Neither", "Maternal RSV vaccine")
)
nir_leaves <- build_period_stack(
  nir_full_year, nir_in_season,
  x_title = "Days since matched index date",
  arms = c("Neither", "Nirsevimab")
)

figure1_km <- (
  mat_leaves$top_curve    + nir_leaves$top_curve    +
    mat_leaves$top_risk     + nir_leaves$top_risk     +
    mat_leaves$bottom_curve + nir_leaves$bottom_curve +
    mat_leaves$bottom_risk  + nir_leaves$bottom_risk
) +
  plot_layout(ncol = 2, heights = c(4, 1, 4, 1))

print(figure1_km)
ggsave("figure1_km_full_year_and_rsv_season_V2_1.pdf", figure1_km, width = 12, height = 12, units = "in", device = "pdf")

build_single_period_figure <- function(fit_obj, x_title, arms, panel_title,
                                       y_max, y_breaks) {
  curve <- extract_km_curve(fit_obj$fit, fit_obj$data, arms = arms, xmax = FOLLOWUP_XMAX)
  p <- plot_final_km(curve, panel_title, x_title, FOLLOWUP_XMAX, arms,
                     y_max = y_max, y_breaks = y_breaks)
  risk <- km_risk_table(fit_obj$data)
  plot_km_with_risk(p, risk, arms, x_title, FOLLOWUP_XMAX)
}

mat_oos_curve <- extract_km_curve(
  maternal_out_of_season$fit, maternal_out_of_season$data,
  arms = c("Neither", "Maternal RSV vaccine"), xmax = FOLLOWUP_XMAX
)
nir_oos_curve <- extract_km_curve(
  nir_out_of_season$fit, nir_out_of_season$data,
  arms = c("Neither", "Nirsevimab"), xmax = FOLLOWUP_XMAX
)
neg_sc <- auto_y_scale(mat_oos_curve, nir_oos_curve)

p_neg_mat <- build_single_period_figure(
  maternal_out_of_season, "Days since birth",
  c("Neither", "Maternal RSV vaccine"),
  panel_title = "Out-of-season (negative control)",
  y_max = neg_sc$y_max, y_breaks = neg_sc$y_breaks
)
p_neg_nir <- build_single_period_figure(
  nir_out_of_season, "Days since matched index date",
  c("Neither", "Nirsevimab"),
  panel_title = "Out-of-season (negative control)",
  y_max = neg_sc$y_max, y_breaks = neg_sc$y_breaks
)

figureS_km_negative_control <- p_neg_mat | p_neg_nir
print(figureS_km_negative_control)
ggsave("figureS_km_out_of_season_negative_control_V2_1.pdf", figureS_km_negative_control,
       width = 12, height = 6, units = "in", device = "pdf")
