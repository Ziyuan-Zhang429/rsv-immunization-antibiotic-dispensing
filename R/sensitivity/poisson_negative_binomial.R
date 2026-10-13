

suppressPackageStartupMessages({
 library(dplyr); library(tidyr); library(purrr); library(lubridate)
 library(splines); library(MASS); library(sandwich); library(lmtest)
 library(ggplot2)
})
select <- dplyr::select
filter <- dplyr::filter
stopifnot(exists("model_data_primary"), exists("model_data_fullyear"),
          exists("model_data_out_of_season"), exists("cox_primary"),
          exists("COVARIATE_TERMS"), exists("FOLLOWUP_MAX"), exists("fit_ag"))
NBP_AGE_KNOTS <- c(7,30,90)
NBP_AGE_BOUNDARY <- c(0,FOLLOWUP_MAX)
NBP_AGE_CUTS <- c(1,3,7,14,21,30,45,60,90,120,150)
NBP_OUTPUT_DIR <- "NB_Poisson_all_periods_results"
NBP_RUN_NB <- TRUE
NBP_RUN_POISSON <- TRUE

run_period <- function(input_data, period, period_slug) {
raw <- input_data %>%
  mutate(.source_row = row_number())
required <- c("infant_study_id", "mom_study_id", "index_date", "tstart",
              "tstop", "event", "exposure_arm", "calendar_month")
stopifnot(all(required %in% names(raw)), nrow(raw) > 0L,
          all(raw$tstop > raw$tstart), all(raw$event %in% c(0, 1)),
          !anyNA(raw$mom_study_id), !anyNA(raw$exposure_arm))
raw <- raw %>% mutate(index_date = as.Date(index_date))
stopifnot(!anyNA(raw$index_date))

split_one <- function(i) {
  row <- raw[i, , drop = FALSE]
  lo <- row$tstart[[1]]
  hi <- row$tstop[[1]]
  birth <- row$index_date[[1]]
  start_date <- birth + floor(lo)
  end_date <- birth + ceiling(hi)
  first_month <- floor_date(start_date, "month")
  last_month <- floor_date(end_date, "month") + months(1)
  month_starts <- seq.Date(as.Date(first_month), as.Date(last_month), by = "month")
  month_cuts <- as.numeric(month_starts - birth) - 1
  cuts <- sort(unique(c(lo, hi, NBP_AGE_CUTS, month_cuts)))
  cuts <- cuts[cuts >= lo & cuts <= hi]
  cuts <- sort(unique(c(lo, cuts, hi)))
  out <- row[rep(1L, length(cuts) - 1L), , drop = FALSE]
  out$tstart <- head(cuts, -1)
  out$tstop <- tail(cuts, -1)
  out$event <- 0L
  out$event[nrow(out)] <- row$event[[1]]
  out$calendar_month <- factor(format(birth + out$tstop, "%Y-%m"))
  out
}
split_dat <- dplyr::bind_rows(lapply(seq_len(nrow(raw)), split_one)) %>%
  mutate(person_days = tstop - tstart,
         age_mid = (tstart + tstop) / 2)

by_source <- split_dat %>%
  group_by(.source_row) %>%
  summarise(events = sum(event), days = sum(person_days), .groups = "drop")
stopifnot(nrow(by_source) == nrow(raw),
          identical(as.numeric(by_source$events), as.numeric(raw$event)),
          isTRUE(all.equal(by_source$days, raw$tstop - raw$tstart,
                           tolerance = 1e-9)),
          all(split_dat$person_days > 0),
          sum(split_dat$event) == sum(raw$event),
          isTRUE(all.equal(sum(split_dat$person_days),
                           sum(raw$tstop - raw$tstart), tolerance = 1e-8)))

split_dat <- split_dat %>%
  mutate(calendar_month = factor(as.character(calendar_month)),
         exposure_arm = factor(as.character(exposure_arm),
                               levels = levels(raw$exposure_arm)))
stopifnot("d_neither" %in% levels(split_dat$exposure_arm))
split_dat$exposure_arm <- relevel(split_dat$exposure_arm, ref = "d_neither")

split_dat <- as.data.frame(split_dat)
rownames(split_dat) <- as.character(seq_len(nrow(split_dat)))

nbp_formula <- as.formula(paste0(
  "event ~ exposure_arm + ns(age_mid, knots = c(",
  paste(NBP_AGE_KNOTS, collapse = ", "),
  "), Boundary.knots = c(", paste(NBP_AGE_BOUNDARY, collapse = ", "),
  ")) + calendar_month + ", COVARIATE_TERMS,
  " + offset(log(person_days))"
))
print(nbp_formula)

fit_poisson <- NULL
fit_nb <- NULL
if (NBP_RUN_POISSON) {
  fit_poisson <- glm(nbp_formula, data = split_dat, family = poisson(link = "log"),
                     model = TRUE, x = TRUE, y = TRUE)
}

nb_fit_method <- NA_character_
nb_theta_bounds <- c(0.05, 1000)
fit_nb_stable <- function(formula, dat, poisson_fit) {
  ordinary <- tryCatch(
    suppressWarnings(MASS::glm.nb(formula, data = dat,
      init.theta = 1.2, control = glm.control(maxit = 100),
      model = TRUE, x = TRUE, y = TRUE)),
    error = function(e) {invisible(NULL); NULL}
  )
  if (!is.null(ordinary) && isTRUE(ordinary$converged) &&
      is.finite(ordinary$theta) && ordinary$theta > 0 &&
      all(is.finite(coef(ordinary)))) {
    return(list(fit = ordinary, method = "glm.nb alternating theta"))
  }
  profile <- function(logtheta) {
    theta <- exp(logtheta)
    fit <- tryCatch(suppressWarnings(glm(formula, data = dat,
      family = MASS::negative.binomial(theta = theta),
      control = glm.control(maxit = 100), model = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged) ||
        any(!is.finite(fitted(fit)))) return(1e100)
    -sum(dnbinom(fit$y, size = theta, mu = fitted(fit), log = TRUE))
  }
  opt <- optimize(profile, interval = log(nb_theta_bounds), tol = 0.01)
  if (!is.finite(opt$objective) || opt$objective >= 1e99)
    stop("NB theta profile failed; inspect model design and interval counts")
  theta <- exp(opt$minimum)
  if (theta < nb_theta_bounds[1] * 1.05 || theta > nb_theta_bounds[2] / 1.05)
    warning("NB theta optimum near search boundary: ", signif(theta, 4),
            "; investigate dispersion and interval-splitting sensitivity")
  fit <- glm(formula, data = dat,
    family = MASS::negative.binomial(theta = theta),
    control = glm.control(maxit = 100), model = TRUE, x = TRUE, y = TRUE)
  if (!isTRUE(fit$converged)) stop("Profiled NB GLM did not converge")
  fit$theta <- theta
  list(fit = fit, method = "bounded theta profile")
}
if (NBP_RUN_NB) {
  nb_result <- fit_nb_stable(nbp_formula, split_dat, fit_poisson)
  fit_nb <- nb_result$fit
  nb_fit_method <- nb_result$method
}

robust_tidy <- function(fit, model_label) {

  used <- match(rownames(model.frame(fit)), rownames(split_dat))
  if (anyNA(used) || length(used) != stats::nobs(fit)) {
    stop("Cannot align fitted rows with mother IDs; check model-frame row names")
  }
  cl <- split_dat$mom_study_id[used]
  stopifnot(length(cl) == stats::nobs(fit), !anyNA(cl))
  vc <- sandwich::vcovCL(fit, cluster = cl, type = "HC1")
  est <- coef(fit)
  se <- sqrt(diag(vc))
  terms <- grep("^exposure_arm", names(est), value = TRUE)
  if (length(terms) != 2L) {
    warning("Expected two exposure contrasts; found: ", paste(terms, collapse = ", "))
  }
  tibble(model = model_label, exposure = terms,
         log_ratio = unname(est[terms]), robust_se = unname(se[terms])) %>%
    mutate(ratio = exp(log_ratio),
           CI_low = exp(log_ratio - qnorm(.975) * robust_se),
           CI_high = exp(log_ratio + qnorm(.975) * robust_se),
           VE_pct = 100 * (1 - ratio),
           VE_CI_low = 100 * (1 - CI_high),
           VE_CI_high = 100 * (1 - CI_low),
           p_value = 2 * pnorm(-abs(log_ratio / robust_se)))
}

count_results <- bind_rows(
  if (!is.null(fit_poisson)) robust_tidy(fit_poisson, "Robust Poisson"),
  if (!is.null(fit_nb)) robust_tidy(fit_nb, "Negative binomial")
)

ag_fit <- if (period == "RSV seasons") cox_primary$all_cause else fit_ag(input_data)
ag_coefs <- coef(ag_fit)
ag_se <- sqrt(diag(vcov(ag_fit)))
ag_terms <- grep("^exposure_arm", names(ag_coefs), value = TRUE)
ag_results <- tibble(model = "Andersen-Gill Cox", exposure = ag_terms,
                     log_ratio = unname(ag_coefs[ag_terms]),
                     robust_se = unname(ag_se[ag_terms])) %>%
  mutate(ratio = exp(log_ratio),
         CI_low = exp(log_ratio - qnorm(.975) * robust_se),
         CI_high = exp(log_ratio + qnorm(.975) * robust_se),
         VE_pct = 100 * (1 - ratio),
         VE_CI_low = 100 * (1 - CI_high),
         VE_CI_high = 100 * (1 - CI_low),
         p_value = 2 * pnorm(-abs(log_ratio / robust_se)))

comparison <- bind_rows(ag_results, count_results) %>%
  mutate(period = period) %>%
  mutate(exposure_label = case_when(
    grepl("maternal", exposure) ~ "Maternal RSVpreF",
    grepl("nirsevimab", exposure) ~ "Nirsevimab",
    TRUE ~ exposure
  ),
  ratio_type = if_else(model == "Andersen-Gill Cox", "HR", "IRR"),
  formatted = sprintf("%.1f%% (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high)) %>%
  select(model, exposure_label, ratio_type, ratio, CI_low, CI_high,
         VE_pct, VE_CI_low, VE_CI_high, p_value, formatted, everything())

pearson_dispersion <- if (!is.null(fit_poisson)) {
  sum(residuals(fit_poisson, type = "pearson")^2) / df.residual(fit_poisson)
} else {
  NA_real_
}

summary_checks <- tibble(
  metric = c("Infants", "Mothers", "Original AG rows", "Split count-model rows",
             "Original antibiotic events", "Split antibiotic events",
             "Original person-days", "Split person-days",
             "Number of calendar months", "Poisson Pearson dispersion",
             "NB theta", "NB fitting method"),
  value = c(n_distinct(raw$infant_study_id), n_distinct(raw$mom_study_id),
            nrow(raw), nrow(split_dat), sum(raw$event), sum(split_dat$event),
            sum(raw$tstop - raw$tstart), sum(split_dat$person_days),
            n_distinct(split_dat$calendar_month), pearson_dispersion,
            if (!is.null(fit_nb)) fit_nb$theta else NA_real_,
            NA_real_)
)

summary_checks <- summary_checks %>% mutate(period = period, .before = 1)
summary_checks <- summary_checks %>% mutate(detail = if_else(metric == "NB fitting method", nb_fit_method, NA_character_))
print(comparison %>% select(exposure_label, model, ratio_type, ratio,
                            CI_low, CI_high, formatted), n = Inf)
print(summary_checks, n = Inf)

if (!is.na(pearson_dispersion) && pearson_dispersion > 2)
period_dir <- file.path(NBP_OUTPUT_DIR, period_slug)
dir.create(period_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(comparison, file.path(period_dir, "AG_NB_Poisson_comparison.csv"),
          row.names = FALSE)
write.csv(summary_checks, file.path(period_dir, "data_checks.csv"),
          row.names = FALSE)

plot_data <- comparison %>%
  mutate(model = factor(model,
                        levels = c("Andersen-Gill Cox", "Robust Poisson",
                                   "Negative binomial")))
p <- ggplot(plot_data, aes(x = VE_pct, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = VE_CI_low, xmax = VE_CI_high), height = 0.16) +
  geom_point(size = 2.4) +
  facet_wrap(~ exposure_label, ncol = 1) +
  labs(x = "Estimated reduction in antibiotic event rate / hazard (%)",
       y = NULL,
       subtitle = "Cox: hazard ratio; count models: incidence rate ratio; 95% cluster-robust CI") +
  theme_bw(base_size = 11)
ggsave(file.path(period_dir, "AG_NB_Poisson_comparison.png"), p,
       width = 9, height = 5.5, dpi = 300)

nbp_analysis_data <- split_dat
nbp_poisson_fit <- fit_poisson
nbp_nb_fit <- fit_nb
nbp_comparison <- comparison
nbp_diagnostics <- summary_checks
list(comparison = comparison, diagnostics = summary_checks)

}

period_inputs <- list(
  list(data = model_data_fullyear$all_cause, label = "Full year", slug = "full_year"),
  list(data = model_data_primary$all_cause, label = "RSV seasons", slug = "rsv_seasons"),
  list(data = model_data_out_of_season$all_cause, label = "Non-RSV seasons", slug = "non_rsv_seasons")
)
period_results <- lapply(period_inputs, function(x) {
  run_period(x$data, x$label, x$slug)
})
all_comparison <- dplyr::bind_rows(lapply(period_results, `[[`, "comparison"))
all_diagnostics <- dplyr::bind_rows(lapply(period_results, `[[`, "diagnostics"))

period_order <- c("Full year", "RSV seasons", "Non-RSV seasons")
model_order <- c("Andersen-Gill Cox", "Robust Poisson", "Negative binomial")
supplement_table <- all_comparison %>%
  mutate(period = factor(period, levels = period_order),
         model = factor(model, levels = model_order)) %>%
  select(period, model, exposure_label, formatted) %>%
  tidyr::pivot_wider(names_from = exposure_label, values_from = formatted) %>%
  arrange(period, model) %>%
  mutate(period = as.character(period), model = as.character(model))
stopifnot(nrow(supplement_table) == 9L,
          all(c("Maternal RSVpreF", "Nirsevimab") %in% names(supplement_table)))

dir.create(NBP_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(supplement_table, file.path(NBP_OUTPUT_DIR, "Supplementary_table_all_periods.csv"), row.names = FALSE)
write.csv(all_comparison, file.path(NBP_OUTPUT_DIR, "All_periods_full_estimates.csv"), row.names = FALSE)
write.csv(all_diagnostics, file.path(NBP_OUTPUT_DIR, "All_periods_data_checks.csv"), row.names = FALSE)
print(supplement_table, n = Inf, width = Inf)
