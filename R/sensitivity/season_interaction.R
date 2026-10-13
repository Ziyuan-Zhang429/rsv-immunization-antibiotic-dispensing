

md_interaction <- clean_and_restrict(
  cp_primary$all_cause,
  restrict_to_season = FALSE,
  keep_out_of_season_only = FALSE
) %>%
  mutate(
    maternal = as.integer(exposure_arm == "a_maternal_vax_only"),
    nirsevimab = as.integer(exposure_arm == "b_nirsevimab_only"),
    season_numeric = as.integer(in_season),
    maternal_x_season = maternal * season_numeric,
    nirsevimab_x_season = nirsevimab * season_numeric
  )

md_interaction$utd_status <- droplevels(md_interaction$utd_status)

stopifnot(
  all(!is.na(md_interaction$in_season)),
  all(as.character(md_interaction$exposure_arm) %in%
        c("d_neither", "a_maternal_vax_only", "b_nirsevimab_only")),
  all(md_interaction$tstop > md_interaction$tstart),
  all(c(0L, 1L) %in% unique(md_interaction$season_numeric))
)

interaction_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~",
  "maternal + nirsevimab + maternal_x_season + nirsevimab_x_season +",
  "strata(calendar_month) +", COVARIATE_TERMS
))

cox_season_interaction <- fit_ag(md_interaction, formula = interaction_formula)

interaction_beta <- coef(cox_season_interaction)
interaction_vcov <- vcov(cox_season_interaction)
interaction_terms_required <- c(
  "maternal",
  "nirsevimab",
  "maternal_x_season",
  "nirsevimab_x_season"
)

stopifnot(
  all(interaction_terms_required %in% names(interaction_beta)),
  all(is.finite(interaction_beta[interaction_terms_required])),
  all(is.finite(
    vcov(cox_season_interaction)[
      interaction_terms_required,
      interaction_terms_required
    ]
  ))
)

linear_contrast <- function(weights, label) {
  w <- setNames(rep(0, length(interaction_beta)), names(interaction_beta))
  w[names(weights)] <- weights
  log_hr <- sum(w * interaction_beta)
  se <- sqrt(drop(t(w) %*% interaction_vcov %*% w))
  z <- log_hr / se
  tibble(
    contrast = label,
    HR = exp(log_hr),
    CI_low = exp(log_hr - qnorm(0.975) * se),
    CI_high = exp(log_hr + qnorm(0.975) * se),
    p_value = 2 * pnorm(-abs(z))
  )
}

season_interaction_table <- bind_rows(
  linear_contrast(c(maternal = 1), "Maternal: non-RSV season"),
  linear_contrast(c(maternal = 1, maternal_x_season = 1), "Maternal: RSV season"),
  linear_contrast(c(maternal_x_season = 1), "Maternal: HR ratio (RSV/non-RSV)"),
  linear_contrast(c(nirsevimab = 1), "Nirsevimab: non-RSV season"),
  linear_contrast(c(nirsevimab = 1, nirsevimab_x_season = 1), "Nirsevimab: RSV season"),
  linear_contrast(c(nirsevimab_x_season = 1), "Nirsevimab: HR ratio (RSV/non-RSV)")
) %>%
  mutate(
    VE_pct = if_else(grepl("HR ratio", contrast), NA_real_, 100 * (1 - HR)),
    VE_CI_low = if_else(grepl("HR ratio", contrast), NA_real_, 100 * (1 - CI_high)),
    VE_CI_high = if_else(grepl("HR ratio", contrast), NA_real_, 100 * (1 - CI_low))
  )
print(season_interaction_table, n = Inf)
print(md_interaction %>% group_by(in_season) %>%
        summarise(person_days = sum(tstop - tstart),
                  events = sum(event),
                  infants_contributing = n_distinct(infant_study_id),
                  .groups = "drop"))

interaction_terms <- c("maternal_x_season", "nirsevimab_x_season")
b <- interaction_beta[interaction_terms]
V <- interaction_vcov[interaction_terms, interaction_terms, drop = FALSE]
wald_chisq <- drop(t(b) %*% solve(V, b))
omnibus_interaction <- tibble(
  test = "Both immunization x season terms = 0",
  chi_square = wald_chisq,
  df = length(interaction_terms),
  p_value = pchisq(wald_chisq, df = length(interaction_terms), lower.tail = FALSE)
)
print(omnibus_interaction)

interaction_formula_flexible <- update(
  interaction_formula_flexible,
  . ~ . - strata(calendar_month) +
    strata(calendar_month, season_numeric)
)

cox_season_interaction_flexible <- fit_ag(
  md_interaction,
  formula = interaction_formula_flexible
)

summary(cox_season_interaction_flexible)
