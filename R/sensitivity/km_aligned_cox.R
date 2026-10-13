

library(splines)

WORK_DIR <- getwd()
SCRIPT_DIR <- getwd()
OUT_DIR <- file.path(WORK_DIR, "reviewer2_comment2_outputs")
RUN_WEIGHTED_NIR <- TRUE
RUN_MATCH_STRATIFIED_NIR <- FALSE
PERIODS <- c("full_year", "in_season", "out_of_season")
MATERNAL_KM_START <- as.Date("2023-11-15")
set.seed(20260809)

AG_FILE <- file.path(SCRIPT_DIR, "R/primary/cox_primary.R")
COUNT_FILE <- file.path(SCRIPT_DIR, "R/primary/exposure_counts.R")
KM_FILE <- file.path(SCRIPT_DIR, "R/figures_tables/km_curves.R")
stopifnot(all(file.exists(c(AG_FILE, COUNT_FILE, KM_FILE))))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
setwd(WORK_DIR)

source(AG_FILE, local = FALSE)
source(COUNT_FILE, local = FALSE)
source(KM_FILE, local = FALSE)

needed <- c("cp_all_cause_recurrent", "nir_base", "cohort_primary",
            "model_data_primary", "model_data_fullyear", "model_data_out_of_season",
            "clean_and_restrict", "COVARIATE_TERMS", "fit_ag", "tidy_ve")
stopifnot(all(vapply(needed, exists, logical(1), inherits = TRUE)))
stopifnot(identical(as.Date(POSITIVITY_START), MATERNAL_KM_START))

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(broom)
library(lubridate)

period_subset <- function(x, period, season_column = "in_season") {
  if (period == "full_year") return(x)
  if (period == "in_season") return(x[!is.na(x[[season_column]]) & x[[season_column]], , drop = FALSE])
  if (period == "out_of_season") return(x[!is.na(x[[season_column]]) & !x[[season_column]], , drop = FALSE])
  stop("Unknown period: ", period)
}

fit_two_arm_ag <- function(md, exposure_col, exposed_level, period,
                           analysis, use_weights = FALSE, match_strata = FALSE) {

  md <- md %>% filter(tstop > tstart, !is.na(.data[[exposure_col]]))
  md[[exposure_col]] <- factor(as.character(md[[exposure_col]]),
                               levels = c("Neither", exposed_level))
  predictors <- all.vars(as.formula(paste("~", COVARIATE_TERMS)))
  missing_cols <- setdiff(c(predictors, "mom_study_id", "calendar_month"), names(md))
  if (length(missing_cols)) stop(analysis, ": missing columns: ", paste(missing_cols, collapse=", "))
  md <- md[complete.cases(md[, predictors, drop = FALSE]), , drop = FALSE]
  md <- md %>% mutate(calendar_month = factor(calendar_month))
  if (use_weights) {
    stopifnot("sw_trunc" %in% names(md), all(is.finite(md$sw_trunc)), all(md$sw_trunc > 0))
  }
  n_events <- sum(md$event)
  n_arms <- n_distinct(md[[exposure_col]])
  n_infants <- n_distinct(md$infant_study_id)
  counts <- md %>% group_by(.data[[exposure_col]]) %>% summarise(
    episodes = if ("episode_id" %in% names(md)) n_distinct(episode_id) else NA_integer_,
    infants = n_distinct(infant_study_id), person_days = sum(tstop - tstart),
    events = sum(event), .groups = "drop") %>%
    rename(arm = all_of(exposure_col)) %>% mutate(analysis = analysis, period = period)
  write.csv(counts, file.path(OUT_DIR, paste0(analysis,"_",period,"_arm_counts.csv")), row.names = FALSE)
  if (n_arms != 2L || n_events < 10) {
    warning("Skipping ", analysis, "/",period,": insufficient arms or events")
    return(tibble(analysis=analysis, period=period, weighted=use_weights,
                  n_infants=n_infants, n_events=n_events, VE_pct=NA_real_,
                  VE_CI_low=NA_real_, VE_CI_high=NA_real_, note="not fit"))
  }
  rhs <- paste0(exposure_col, " + strata(calendar_month) + ", COVARIATE_TERMS)
  if (match_strata) rhs <- paste(rhs, "+ strata(match_id)")
  f <- as.formula(paste("Surv(tstart,tstop,event) ~",rhs))
  args <- list(formula=f, data=md, ties="efron", cluster=md$mom_study_id,
               model=TRUE, control=coxph.control(iter.max=35))
  if (use_weights) args$weights <- md$sw_trunc
  fit <- do.call(coxph, args)
  term <- paste0(exposure_col, exposed_level)
  tidy <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>% filter(.data$term == !!term)
  if (nrow(tidy) != 1L) stop("Exposure coefficient missing in ", analysis, " / ", period)
  saveRDS(fit, file.path(OUT_DIR, paste0(analysis,"_",period,"_fit.rds")))
  tibble(analysis=analysis, period=period, weighted=use_weights,
         n_infants=n_infants,
         n_episodes=if ("episode_id" %in% names(md)) n_distinct(md$episode_id) else NA_integer_,
         n_events=n_events, HR=tidy$estimate, HR_low=tidy$conf.low, HR_high=tidy$conf.high,
         VE_pct=100*(1-tidy$estimate), VE_CI_low=100*(1-tidy$conf.high),
         VE_CI_high=100*(1-tidy$conf.low), p_value=tidy$p.value,
         note=if (match_strata) "matched sets stratified" else "calendar-month stratified")
}

reference <- list(in_season=model_data_primary$all_cause,
                  full_year=model_data_fullyear$all_cause,
                  out_of_season=model_data_out_of_season$all_cause)
orig <- imap_dfr(reference, function(md, period) {
  fit <- fit_ag(md)
  broom::tidy(fit, conf.int=TRUE, exponentiate=TRUE) %>%
    filter(grepl("^exposure_arm",term), !grepl(":",term)) %>%
    transmute(analysis="Original_AG", period=period,
              intervention=case_when(grepl("maternal",term) ~ "Maternal RSVpreF",
                                     grepl("nirsevimab",term) ~ "Nirsevimab", TRUE ~ term),
              n_infants=n_distinct(md$infant_study_id),n_events=sum(md$event),
              HR=estimate,HR_low=conf.low,HR_high=conf.high,
              VE_pct=100*(1-estimate),VE_CI_low=100*(1-conf.high),
              VE_CI_high=100*(1-conf.low))
})
write.csv(orig, file.path(OUT_DIR,"original_AG_reference.csv"),row.names=FALSE)

restricted_three_arm <- map_dfr(PERIODS, function(period) {
  md <- reference[[period]] %>% filter(index_date >= MATERNAL_KM_START)
  fit <- fit_ag(md)
  broom::tidy(fit, conf.int=TRUE, exponentiate=TRUE) %>%
    filter(grepl("^exposure_arm", term), !grepl(":",term)) %>%
    transmute(analysis="M0_birth_Nov15_three_arm",period=period,
      intervention=case_when(grepl("maternal",term) ~ "Maternal RSVpreF",
                             grepl("nirsevimab",term) ~ "Nirsevimab", TRUE ~ term),
      n_infants=n_distinct(md$infant_study_id),n_events=sum(md$event),
      HR=estimate,HR_low=conf.low,HR_high=conf.high,
      VE_pct=100*(1-estimate),VE_CI_low=100*(1-conf.high),
      VE_CI_high=100*(1-conf.low))
})
write.csv(restricted_three_arm,file.path(OUT_DIR,"M0_birth_restricted_three_arm.csv"),row.names=FALSE)

mat_cp <- cp_all_cause_recurrent %>%
  filter(index_date >= MATERNAL_KM_START) %>%
  left_join(cohort_primary %>% transmute(infant_study_id,
              subsequent_mab_day = if_else(mab_history, pmax(mab_days,0),Inf)),
            by="infant_study_id") %>%
  mutate(arm2=case_when(exposure_arm=="d_neither" ~ "Neither",
                        exposure_arm=="a_maternal_vax_only" ~ "Maternal RSVpreF",
                        TRUE ~ NA_character_))
stopifnot(!anyNA(mat_cp$subsequent_mab_day))

make_maternal <- function(period, censor_at_mab=FALSE) {
  x <- period_subset(mat_cp, period)
  if (censor_at_mab) {

    x <- x %>%
      filter(tstart < subsequent_mab_day) %>%
      mutate(tstop_old=tstop, tstop=pmin(tstop,subsequent_mab_day),
             event=if_else(tstop_old <= subsequent_mab_day, event,0L),
             event=if_else(is.finite(subsequent_mab_day) &
                              tstop_old == subsequent_mab_day,0L,event)) %>%
      select(-tstop_old)
  }
  x <- x %>% filter(arm2 %in% c("Neither","Maternal RSVpreF"),tstop>tstart)

  apply_unknown_exclusions(x)
}

maternal_results <- map_dfr(PERIODS, function(period) {
  bind_rows(
    fit_two_arm_ag(make_maternal(period,FALSE),"arm2","Maternal RSVpreF",period,"M1_birth_Nov15"),
    fit_two_arm_ag(make_maternal(period,TRUE), "arm2","Maternal RSVpreF",period,"M2_birth_Nov15_censor_mab")
  )
})

nir_episodes <- nir_base %>% distinct(episode_id, infant_study_id, match_id, exposure_label,
                            index_age, raw_max_followup, km_index_date, sw_trunc)
stopifnot(nrow(nir_episodes)==n_distinct(nir_base$episode_id),
          all(nir_base$tstop>nir_base$tstart),
          all(nir_base$tstop<=nir_base$raw_max_followup + 1e-8),
          all(nir_base$event %in% c(0L,1L)),
          !anyNA(nir_base$sw_trunc))

expand_age_30 <- function(df) {
  vaccine_cols <- names(cohort_primary)[grepl("^(HBV_|rota_|DTaP_|Hib_|PCV_|IPV_)",names(cohort_primary))]
  stopifnot(length(vaccine_cols)>0)
  vax <- cohort_primary %>% select(infant_study_id,all_of(vaccine_cols))
  x <- df %>% mutate(.row_id=row_number(),
                      age_start=index_age+tstart, age_stop=index_age+tstop,
                      .cuts=map2(age_start,age_stop,function(a,b) {
                        inner <- seq(30,150,30)
                        sort(unique(c(a,inner[inner>a & inner<b],b)))
                      }))
  x <- x %>% unnest_longer(.cuts,values_to="cut") %>%
    group_by(.row_id) %>% arrange(cut,.by_group=TRUE) %>%
    mutate(new_tstart=cut,new_tstop=lead(cut)) %>% ungroup() %>%
    filter(!is.na(new_tstop),new_tstop>new_tstart) %>%
    mutate(event=if_else(abs(new_tstop-age_stop)<1e-8,event,0L),
           tstart=new_tstart,tstop=new_tstop) %>%
    select(-cut,-new_tstart,-new_tstop,-age_start,-age_stop,-.row_id,-utd_status) %>%
    left_join(vax,by="infant_study_id") %>%
    mutate(interval_start=tstart) %>% build_utd_at(interval_start) %>%
    select(-all_of(vaccine_cols),-interval_start)
  x
}

nir_ag <- expand_age_30(nir_base) %>%
  mutate(arm2=as.character(exposure_label),

         calendar_date_start=km_index_date + (tstart-index_age) + 1,
         calendar_month=floor_date(calendar_date_start,"month"))
stopifnot(all(nir_ag$tstart >= nir_ag$index_age),
          all(nir_ag$tstop <= 180+1e-8),
          sum(nir_ag$event)==sum(nir_base$event),
          all(nir_ag$arm2 %in% c("Neither","Nirsevimab")),
          all(nir_ag$tstop>nir_ag$tstart))

diagnostics <- tibble(
  metric=c("KM matched episodes", "Unique infants", "Unique treated matches",
           "Control episodes", "Reused control infants", "Event rows KM", "Event rows AG"),
  value=c(n_distinct(nir_base$episode_id),n_distinct(nir_base$infant_study_id),
          n_distinct(nir_base$match_id),
          n_distinct(nir_base$episode_id[nir_base$exposure_label=="Neither"]),
          nir_base %>% filter(exposure_label=="Neither") %>%
            distinct(episode_id,infant_study_id) %>% count(infant_study_id) %>%
            summarise(n=sum(n>1)) %>% pull(n),
          sum(nir_base$event),sum(nir_ag$event)))
write.csv(diagnostics,file.path(OUT_DIR,"nirsevimab_matched_diagnostics.csv"),row.names=FALSE)

nir_results <- map_dfr(PERIODS, function(period) {
  dat <- period_subset(nir_ag,period,season_column="seg_in_season")

  results <- fit_two_arm_ag(dat,"arm2","Nirsevimab",period,"N1_matched_postindex")
  if (RUN_WEIGHTED_NIR) {
    results <- bind_rows(results,fit_two_arm_ag(dat,"arm2","Nirsevimab",period,
                         "N2_matched_postindex_IPTW",use_weights=TRUE))
  }
  if (RUN_MATCH_STRATIFIED_NIR) {
    results <- bind_rows(results,fit_two_arm_ag(dat,"arm2","Nirsevimab",period,
                         "N3_match_stratified",match_strata=TRUE))
  }
  results
})

all_results <- bind_rows(maternal_results,nir_results) %>%
  mutate(intervention=if_else(grepl("^M",analysis),"Maternal RSVpreF","Nirsevimab"))
write.csv(all_results,file.path(OUT_DIR,"KM_aligned_AG_sensitivity_estimates.csv"),row.names=FALSE)

compare <- all_results %>%
  left_join(orig %>% select(period,intervention,primary_VE=VE_pct,
                            primary_CI_low=VE_CI_low,primary_CI_high=VE_CI_high),
            by=c("period","intervention")) %>%
  mutate(VE_difference_pp=VE_pct-primary_VE) %>%
  arrange(intervention,period,analysis)
write.csv(compare,file.path(OUT_DIR,"KM_aligned_vs_primary_comparison.csv"),row.names=FALSE)
print(compare,n=Inf,width=Inf)

writeLines(c(
  "Reviewer 2 Comment 2(d) analyses -- methodological notes",
  "M0: original three-arm AG with birth date >= 2023-11-15; separates eligibility from contrast changes.",
  "M1: birth date >= 2023-11-15; source AG counting process; restrict to two-arm maternal vs neither person-time.",
  "M2: M1 plus censor at nirsevimab dose for previously unexposed infants, including removal of dose-day antibiotic events.",
  "N1: exactly the supplied KM's matched nirsevimab episodes, control switching/censoring, and post-index event process; recurrent events; delayed-entry AGE Cox with calendar month strata and AG covariates.",
  "N2: N1 with the supplied KM's episode-level stabilized truncated IPTW.",
  "N1/N2 include repeated copies of reused control infants; cluster by mother (not episode), but statistical interpretation requires care.",
  "N2 includes both IPT weighting and covariate adjustment; it is a sensitivity specification, not identical to the original KM estimator.",
  "KM nirsevimab original uses days since matched index; sensitivity here uses attained age with delayed entry.",
  "Seasonal analysis filters calendar-season segments before AG estimation; recurrent events in retained segments are counted.",
  "The maternal KM weights are row-level, fitted after first-event restriction; hence not attached naively to recurrent AG rows.",
  "These analyses are not strictly identical estimands; report counts, direction, precision, and changes alongside primary AG.",
  "No SAS data were supplied to this chat, so these output results cannot be validated here."
),file.path(OUT_DIR,"README_interpretation.txt"))
