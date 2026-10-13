

library(dplyr)
library(splines)
library(glmnet)
library(survival)
library(tidyr)
library(purrr)

stopifnot(all(vapply(
  c("cohort_primary", "cp_all_cause_recurrent", "COVARIATE_TERMS",
    "fit_two_arm_ag", "apply_unknown_exclusions", "OUT_DIR", "PERIODS"),
  exists, logical(1), inherits = TRUE
)))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

mat_cp <- cp_all_cause_recurrent %>%
  filter(index_date >= as.Date("2023-11-15")) %>%
  mutate(
    maternal_arm = case_when(
      as.character(exposure_arm) == "a_maternal_vax_only" ~ "Maternal RSVpreF",
      as.character(exposure_arm) == "d_neither" ~ "Neither",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(maternal_arm)) %>%
  apply_unknown_exclusions()

baseline_terms <- gsub(" ?\\+ ?utd_status", "", COVARIATE_TERMS)
baseline_vars <- all.vars(as.formula(paste("~", baseline_terms)))

eligible_ids <- mat_cp %>% distinct(infant_study_id)
baseline <- cohort_primary %>%
  semi_join(eligible_ids, by = "infant_study_id") %>%
  mutate(
    maternal_arm = if_else(rsv_valid, "Maternal RSVpreF", "Neither"),
    birth_day_numeric = as.numeric(index_date)
  )

missing_vars <- setdiff(c("infant_study_id", "index_date",
                          "birth_day_numeric", baseline_vars), names(baseline))
if (length(missing_vars)) {
  stop("Missing baseline PS variables: ", paste(missing_vars, collapse = ", "))
}

baseline <- baseline[complete.cases(
  baseline[, c("birth_day_numeric", baseline_vars), drop = FALSE]
), , drop = FALSE]

stopifnot(!anyDuplicated(baseline$infant_study_id),
          n_distinct(baseline$maternal_arm) == 2L)

den_formula <- as.formula(paste(
  "~ splines::ns(birth_day_numeric, df=3) +", baseline_terms
))
num_formula <- ~ splines::ns(birth_day_numeric, df=3)

matrix_no_intercept <- function(f, data) {
  mm <- model.matrix(f, data = data)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}
Xden <- matrix_no_intercept(den_formula, baseline)
Xnum <- matrix_no_intercept(num_formula, baseline)
y <- as.integer(baseline$maternal_arm == "Maternal RSVpreF")
set.seed(20260809)
den_fit <- cv.glmnet(Xden, y, family="binomial", alpha=0, nfolds=5,
                     standardize=TRUE)
set.seed(20260809)
num_fit <- cv.glmnet(Xnum, y, family="binomial", alpha=0, nfolds=5,
                     standardize=TRUE)
pden <- as.numeric(predict(den_fit, newx=Xden, s="lambda.1se", type="response"))
pnum <- as.numeric(predict(num_fit, newx=Xnum, s="lambda.1se", type="response"))
pden <- pmin(pmax(pden, 0.01), 0.99)
pnum <- pmin(pmax(pnum, 0.01), 0.99)

baseline <- baseline %>%
  mutate(
    ps_den = pden, ps_num = pnum,
    sw_raw = if_else(y == 1L, ps_num / ps_den,
                    (1-ps_num)/(1-ps_den))
  )
bounds <- quantile(baseline$sw_raw, c(0.01, 0.99), na.rm=TRUE)
baseline <- baseline %>%
  mutate(sw_trunc = pmin(pmax(sw_raw, bounds[[1]]), bounds[[2]]))

weighted_cp <- mat_cp %>%
  inner_join(
    baseline %>% select(infant_study_id, sw_trunc),
    by = "infant_study_id", relationship = "many-to-one"
  ) %>%
  mutate(arm2 = maternal_arm)
stopifnot(n_distinct(weighted_cp$infant_study_id) == nrow(baseline),
          all(is.finite(weighted_cp$sw_trunc)),
          all(weighted_cp$sw_trunc > 0))

maternal_iptw_results <- purrr::map_dfr(PERIODS, function(period) {
  dat <- switch(period,
    full_year = weighted_cp,
    in_season = weighted_cp %>% filter(in_season),
    out_of_season = weighted_cp %>% filter(!in_season)
  )
  fit_two_arm_ag(
    md = dat, exposure_col = "arm2",
    exposed_level = "Maternal RSVpreF",
    period = period,
    analysis = "M3_birth_Nov15_baseline_IPTW",
    use_weights = TRUE
  )
})

write.csv(maternal_iptw_results,
          file.path(OUT_DIR, "M3_maternal_IPTW_AG_estimates.csv"),
          row.names=FALSE)
write.csv(
  baseline %>% group_by(maternal_arm) %>%
    summarise(n_infants=n(), mean_weight=mean(sw_trunc),
              median_weight=median(sw_trunc),
              min_weight=min(sw_trunc), max_weight=max(sw_trunc),
              mean_ps_den=mean(ps_den), .groups="drop"),
  file.path(OUT_DIR, "M3_maternal_weight_diagnostics.csv"),
  row.names=FALSE
)
write.csv(
  tibble(
    n_baseline_infants=nrow(baseline),
    n_recurrent_events=sum(weighted_cp$event),
    n_weighted_rows=nrow(weighted_cp),
    trunc_lower=as.numeric(bounds[[1]]),
    trunc_upper=as.numeric(bounds[[2]]),
    n_floored_den=sum(pden %in% c(0.01,0.99)),
    n_floored_num=sum(pnum %in% c(0.01,0.99))
  ),
  file.path(OUT_DIR, "M3_maternal_weight_overview.csv"),
  row.names=FALSE
)
print(maternal_iptw_results, n=Inf, width=Inf)
