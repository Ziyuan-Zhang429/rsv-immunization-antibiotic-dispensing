

library(survival)
library(dplyr)
library(purrr)
library(tidyr)
library(future)
library(furrr)
library(ggplot2)
library(patchwork)
library(scales)

BOOTSTRAP_REPLICATES <- 200L
n_workers <- min(8, max(1, future::availableCores() - 2))
stopifnot(FOLLOWUP_MAX == 180,
          exists("model_data_fullyear"), exists("model_data_primary"),
          exists("analytic_infants"), exists("arm_counts"))
truncate_cp <- function(md, cutoff) md %>%
  filter(tstart < cutoff) %>%
  mutate(event = if_else(tstop > cutoff, 0L, event), tstop = pmin(tstop, cutoff)) %>%
  filter(tstop > tstart)

fit_ag <- function(md, cutoff, formula = ag_formula) {
  d <- truncate_cp(md, cutoff)
  list(fit = coxph(formula, data = d, cluster = mom_study_id,
                   ties = "efron", model = TRUE), data = d)
}

baseline_funs <- function(fit) suppressWarnings(basehaz(fit, centered = FALSE)) %>%
  mutate(stratum_key = sub("^stratum=", "", as.character(strata))) %>%
  arrange(stratum_key, time) %>%
  split(.$stratum_key) %>%
  map(~stepfun(.x$time, c(0, .x$hazard), right = FALSE))

expected_events <- function(fit, md, H0, arm = NULL, id_col = "infant_study_id") {
  d <- md
  if (!is.null(arm)) d$exposure_arm <- factor(arm, levels = levels(md$exposure_arm))
  lp <- predict(fit, newdata = d, type = "lp", reference = "zero")
  if (any(!is.finite(lp)) || any(!is.finite(exp(lp)))) stop("Non-finite linear predictor.")
  h_start <- h_stop <- numeric(nrow(d))
  s <- as.character(d$calendar_month)
  for (k in names(H0)) {
    i <- which(s == k)
    if (length(i)) { h_start[i] <- H0[[k]](d$tstart[i]); h_stop[i] <- H0[[k]](d$tstop[i]) }
  }
  d$.expected <- (h_stop - h_start) * exp(lp)
  d %>% group_by(.data[[id_col]]) %>%
    summarise(expected_events = sum(.expected), .groups = "drop")
}

standardize_att <- function(fit, md, H0, exposed_arm, label,
                            id_col = "infant_study_id") {
  ids <- md %>% filter(exposure_arm == exposed_arm) %>%
    distinct(.data[[id_col]]) %>% pull(.data[[id_col]])
  if (!length(ids)) stop("No recipients found for ", label, ".")
  d <- md %>% filter(.data[[id_col]] %in% ids)
  obs <- expected_events(fit, d, H0, id_col = id_col) %>%
    rename(expected_protection = expected_events)
  no <- expected_events(fit, d, H0, arm = "d_neither", id_col = id_col) %>%
    rename(expected_no_protection = expected_events)
  z <- left_join(no, obs, by = id_col) %>%
    mutate(diff = expected_no_protection - expected_protection)
  tibble(exposure = label, n_recipients_contributing = length(ids),
         mu_no_protection = mean(z$expected_no_protection),
         mu_protection = mean(z$expected_protection),
         courses_averted_per_child = mean(z$diff),
         courses_averted_per_100 = 100 * mean(z$diff))
}

standardize_both <- function(fit, md, H0, id_col = "infant_study_id") bind_rows(
  standardize_att(fit, md, H0, "a_maternal_vax_only", "Maternal RSV vaccine", id_col),
  standardize_att(fit, md, H0, "b_nirsevimab_only", "Nirsevimab", id_col)
)

run_cutoffs <- function(md) map_dfr(c(90, 180), function(cutoff) {
  x <- fit_ag(md, cutoff); H0 <- baseline_funs(x$fit)
  standardize_both(x$fit, x$data, H0) %>% mutate(cutoff_days = cutoff, .before = 1)
})

bootstrap_courses <- function(md, B = BOOTSTRAP_REPLICATES,
                              id_col = "infant_study_id", cluster_col = "mom_study_id",
                              seed = 20260817) {
  set.seed(seed); moms <- unique(md[[cluster_col]])
  future_map_dfr(seq_len(B), function(b) {
    sampled <- sample(moms, length(moms), replace = TRUE)
    boot_md <- map_dfr(seq_along(sampled), function(i) md %>%
                         filter(.data[[cluster_col]] == sampled[i]) %>%
                         mutate(.boot_infant_id = paste0(.data[[id_col]], "_boot", i),
                                .boot_mom_id = paste0(.data[[cluster_col]], "_boot", i)))
    out <- tryCatch(map_dfr(c(90, 180), function(cutoff) {
      d <- truncate_cp(boot_md, cutoff)
      fit <- coxph(ag_formula, data = d, cluster = .boot_mom_id,
                   ties = "efron", model = TRUE)
      standardize_both(fit, d, baseline_funs(fit), id_col = ".boot_infant_id") %>%
        mutate(cutoff_days = cutoff, .before = 1)
    }), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    mutate(out, rep = b)
  }, .options = furrr_options(seed = TRUE))
}

summarize_bootstrap <- function(point, boot) boot %>%
  group_by(exposure, cutoff_days) %>%
  summarise(ci_low = quantile(courses_averted_per_100, .025, na.rm = TRUE),
            ci_high = quantile(courses_averted_per_100, .975, na.rm = TRUE),
            n_boot = n(), .groups = "drop") %>%
  left_join(point, by = c("exposure", "cutoff_days")) %>%
  arrange(exposure, cutoff_days)

run_analysis <- function(md, panel, seed) {
  point <- run_cutoffs(md); print(point, n = Inf)
  plan(multisession, workers = n_workers)
  boot <- bootstrap_courses(md, seed = seed)
  plan(sequential)
  ci <- summarize_bootstrap(point, boot); print(ci, n = Inf)
  ci %>% transmute(panel, exposure,
                   followup = factor(paste0("1-", cutoff_days, " days"),
                                     levels = c("1-90 days", "1-180 days")),
                   n_recipients_contributing, estimate = courses_averted_per_100, ci_low, ci_high)
}

figure4_data <- bind_rows(
  run_analysis(model_data_fullyear$all_cause, "Full year", 20260817),
  run_analysis(model_data_primary$all_cause, "RSV season", 20260818)
) %>% mutate(panel = factor(panel, levels = c("Full year", "RSV season"))) %>%
  arrange(panel, factor(exposure, levels = c("Maternal RSV vaccine", "Nirsevimab")), followup)

invisible(NULL); print(figure4_data, n = Inf)

make_courses_averted_plot <- function(df) {
  FONT <- 11; PAD <- 2
  panel_order <- c("Full year", "RSV season")
  exposure_order <- c("Maternal RSV vaccine", "Nirsevimab")
  headers <- crossing(panel = panel_order, exposure = exposure_order) %>%
    mutate(followup = NA_character_, n_recipients_contributing = NA_integer_, estimate = NA_real_,
           ci_low = NA_real_, ci_high = NA_real_, row_type = "exposure")
  panel_headers <- tibble(panel = panel_order, exposure = NA_character_, followup = NA_character_,
                          n_recipients_contributing = NA_integer_, estimate = NA_real_, ci_low = NA_real_,
                          ci_high = NA_real_, row_type = "panel")
  dat <- df %>% mutate(panel = as.character(panel), exposure = as.character(exposure),
                       followup = as.character(followup), row_type = "data")
  plot_df <- map_dfr(panel_order, function(p) bind_rows(
    filter(panel_headers, panel == p),
    map_dfr(exposure_order, function(e) bind_rows(
      filter(headers, panel == p, exposure == e),
      filter(dat, panel == p, exposure == e)
    ))
  )) %>%
    mutate(
      row_id = row_number(),
      y_pos = -row_id,
      estimate_annotation = if_else(
        row_type == "data",
        sprintf("%.2f (%.2f, %.2f)", estimate, ci_low, ci_high),
        NA_character_
      )
    )

  n_rows <- nrow(plot_df)
  Y <- c(-(n_rows + .5), -.5)

  shade <- plot_df %>%
    filter(row_id %% 2 == 0)
  transparent <- theme(panel.background = element_blank(),
                       plot.background = element_rect(fill = "transparent", colour = NA))
  LABEL_X <- 0; FOLLOWUP_X <- .18; N_X <- 2.25; TABLE_MIN <- -.25; TABLE_MAX <- 3.35

  p_bg <- ggplot(shade, aes(y = y_pos)) +
    geom_rect(aes(ymin = y_pos - .5, ymax = y_pos + .5, xmin = -Inf, xmax = Inf), fill = "grey88", alpha = .5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) + theme_void() +
    theme(plot.margin = margin(5.5, 0, 5.5, 0), plot.background = element_rect(fill = "white", colour = NA))

  p_table <- ggplot(plot_df, aes(y = y_pos)) +
    geom_text(data = filter(plot_df, row_type == "panel"), aes(x = LABEL_X, label = panel),
              hjust = 0, fontface = "bold", size = FONT / .pt) +
    geom_text(data = filter(plot_df, row_type == "exposure"), aes(x = LABEL_X + FOLLOWUP_X, label = exposure),
              hjust = 0, fontface = "bold", size = FONT / .pt) +
    geom_text(data = filter(plot_df, row_type == "data"), aes(x = LABEL_X + 2 * FOLLOWUP_X, label = followup),
              hjust = 0, size = FONT / .pt) +
    geom_text(data = filter(plot_df, row_type == "data"), aes(x = N_X, label = comma(n_recipients_contributing)),
              hjust = 0, size = FONT / .pt) +
    scale_x_continuous(limits = c(TABLE_MIN, TABLE_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) + theme_void(base_size = FONT) +
    theme(plot.margin = margin(5.5, 8, 5.5, 5.5)) + transparent

  table_header <- ggplot() +
    annotate("text", x = LABEL_X, y = .6, label = "Exposure arm", hjust = 0, fontface = "bold", size = FONT / .pt) +
    annotate("text", x = N_X, y = .6, label = "Recipients\nstandardized", hjust = 0,
             fontface = "bold", size = FONT / .pt) +
    scale_x_continuous(limits = c(TABLE_MIN, TABLE_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 1.5), expand = c(0, 0)) + theme_void() +
    theme(plot.margin = margin(5.5, PAD, 0, 5.5)) + transparent

  p_forest <- ggplot(filter(plot_df, row_type == "data"), aes(x = estimate, y = y_pos)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = .5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = .8, color = "black") +
    geom_point(size = 3.2, shape = 15, color = "black") +
    scale_x_continuous(limits = c(-0.5, 5), breaks = pretty(c(0, 5), n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) +
    labs(x = "Antibiotic courses averted per 100 intervention recipients", y = NULL) +
    theme_minimal(base_size = FONT) +
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = FONT, color = "black"),
          axis.text.x = element_text(size = FONT, color = "black"), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.x = element_line(color = "black", linewidth = .6),
          axis.ticks.x = element_line(color = "black", linewidth = .6), legend.position = "none",
          plot.margin = margin(5.5, 8, 5.5, 8)) + transparent
  forest_header <- plot_spacer() + transparent

  p_est <- ggplot(plot_df, aes(y = y_pos)) +
    geom_text(data = filter(plot_df, row_type == "data"), aes(x = 0, label = estimate_annotation),
              hjust = 0, size = FONT / .pt) +
    scale_x_continuous(limits = c(0, 2.4), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) + theme_void(base_size = FONT) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, PAD)) + transparent
  est_header <- ggplot() +
    annotate("text", x = 0, y = .6, label = "Courses averted\n(95% CI)", hjust = 0,
             fontface = "bold", size = FONT / .pt) +
    scale_x_continuous(limits = c(0, 2.4), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 1.5), expand = c(0, 0)) + theme_void() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 8)) + transparent

  design <- c(area(1,1,1,1), area(1,2,1,2), area(1,3,1,3), area(2,1,2,3),
              area(2,1,2,1), area(2,2,2,2), area(2,3,2,3))
  table_header + forest_header + est_header + p_bg + p_table + p_forest + p_est +
    plot_layout(design = design, heights = c(1, n_rows), widths = c(1.35, .9, .85)) +
    plot_annotation(
      title = "Adjusted Antibiotic Courses Averted per 100 Intervention Recipients",
      
      theme = theme(plot.title = element_text(size = FONT),
                    
                    plot.background = element_rect(fill = "white", colour = NA)))
}

fig4 <- make_courses_averted_plot(figure4_data)
print(fig4)
ggsave("figure4_courses_averted_full_year_rsv_season_0826.pdf", fig4,
       width = 7.5, height = 5.2, units = "in", device = "pdf")
