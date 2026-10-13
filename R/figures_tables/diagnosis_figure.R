

stopifnot(
  exists("arm_counts"), exists("format_events_denom"),
  exists("fig_denom_header"), exists("fig_denom_note"),
  exists("diagnosis_flags"), exists("diagnosis_interpretation"),
  exists("build_diagnosis_events"), exists("clean_and_restrict"),
  exists("build_cp"), exists("COVARIATE_TERMS"),
  exists("empty_ve_row"), exists("tidy_ve"), exists("blocks_primary"),
  "dx_group" %in% names(diagnosis_interpretation)
)

library(tidyverse)
library(survival)
library(broom)
library(ggplot2)
library(patchwork)
library(scales)

add_calendar_3month <- function(df) {
  df %>% mutate(
    calendar_3month = factor(paste0(
      format(calendar_date_start, "%Y"), "_Q",
      (as.integer(format(calendar_date_start, "%m")) - 1L) %/% 3L + 1L
    ))
  )
}

ag_formula_3month <- as.formula(paste(
  "Surv(tstart, tstop, event) ~ exposure_arm + strata(calendar_3month) +", COVARIATE_TERMS
))

restrict_primary_event_days <- function(events_df) {
  stopifnot(exists("cohort_primary"), exists("FOLLOWUP_MAX"))
  events_df %>%
    filter(!is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX) %>%
    semi_join(cohort_primary %>% distinct(infant_study_id), by = "infant_study_id") %>%
    distinct(infant_study_id, abx_days)
}

figure_empty_ve_row <- function(outcome_label, status, n_events_raw = NA_integer_,
                                n_events_model = NA_integer_) {
  empty_ve_row(outcome_label, status, n_events_raw, n_events_model) %>%
    mutate(n_infants = NA_integer_, person_days = NA_real_,
           person_years = NA_real_, n_events = NA_integer_)
}

DIAGNOSIS_LABELS <- c(
  "Bronchitis/bronchiolitis"          = "bronchitis_yn",
  "Fever"                             = "fever_yn",
  "Influenza"                         = "influenza_yn",
  "Bacterial infection (misc.)"       = "bacterial_infection_yn",
  "Other infection (misc.)"           = "other_infection_yn",
  "Acute/suppurative otitis media"    = "sup_om_yn",
  "Non-suppurative otitis media"      = "nsup_om",
  "Pharyngitis"                       = "pharyngitis_yn",
  "Pneumonia"                         = "pneumonia_yn",
  "Sinusitis"                         = "sinusitis_yn",
  "Urinary tract infection"           = "UTI_yn",
  "Skin/cutaneous/mucosal condition"  = "skin_cond_yn",
  "Viral upper respiratory infection" = "URI_yn"
)

DIAGNOSES_TO_DROP <- c("asthma_allergy_yn", "nsup_om", "bacterial_infection_yn",
                       "sinusitis_yn", "UTI_yn", "skin_cond_yn",
                       "other_infection_yn", "fever_yn")

diagnosis_flags_fig <- setdiff(diagnosis_flags, DIAGNOSES_TO_DROP)

fit_single_outcome_panel_supp <- function(outcome_events, outcome_label,
                                     formula = ag_formula_3month,
                                     blocks = blocks_primary,
                                     season_mode = c("in_season", "full_year", "out_season"),
                                     min_events_total = 50,
                                     min_events_both_arm = 5) {

  season_mode <- match.arg(season_mode)

  outcome_events <- restrict_primary_event_days(outcome_events)
  n_ev_raw <- nrow(outcome_events)
  if (n_ev_raw < min_events_total) {
    return(figure_empty_ve_row(outcome_label, "skipped_sparse", n_ev_raw))
  }

  cp <- build_cp(outcome_events, blocks = blocks)

  md <- switch(season_mode,
               in_season  = clean_and_restrict(cp, restrict_to_season = TRUE),
               full_year  = clean_and_restrict(cp, restrict_to_season = FALSE),
               out_season = clean_and_restrict(cp, keep_out_of_season_only = TRUE)
  ) %>% add_calendar_3month()

  has_both <- "c_both" %in% levels(md$exposure_arm)
  n_events_both_arm <- if (has_both) {
    md %>% filter(exposure_arm == "c_both") %>% summarise(n = sum(event)) %>% pull(n)
  } else NA_integer_

  drop_both <- has_both && n_events_both_arm < min_events_both_arm
  if (drop_both) {
    md <- md %>%
      filter(exposure_arm != "c_both") %>%
      mutate(exposure_arm = droplevels(exposure_arm))
  }

  counts <- arm_counts(md, outcome_label) %>%
    select(exposure_arm, exposure, n_infants = n_infants_contributing,
           person_days, person_years, n_events)

  tryCatch({
    fit <- coxph(formula, data = md, cluster = mom_study_id, ties = "efron", model = TRUE)
    ve_rows <- tidy_ve(fit, outcome_label) %>%
      select(-any_of(c("n_infants", "n_events", "person_days", "person_years",
                       "exposure_arm"))) %>%
      left_join(counts %>% select(exposure, n_infants, person_days,
                                  person_years, n_events), by = "exposure") %>%
      mutate(
        n_events_raw   = n_ev_raw,
        n_events_model = sum(md$event),
        fit_status     = if_else(drop_both, "fit_both_dropped_sparse", "fit")
      )

    neither_counts <- counts %>% filter(exposure_arm == "d_neither")
    neither_row <- if (nrow(neither_counts) == 1) {
      tibble(
        outcome        = outcome_label,
        exposure       = "exposure_armd_neither",
        HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
        n_infants      = neither_counts$n_infants,
        person_days    = neither_counts$person_days,
        person_years   = neither_counts$person_years,
        n_events       = neither_counts$n_events,
        n_events_raw   = n_ev_raw,
        n_events_model = sum(md$event),
        fit_status     = if_else(drop_both, "fit_both_dropped_sparse", "fit")
      )
    } else NULL

    bind_rows(ve_rows, neither_row)
  }, error = function(e) {
    figure_empty_ve_row(outcome_label, paste("error:", conditionMessage(e)),
                 n_ev_raw, sum(md$event))
  })
}

diagnosis_panel_results_supp <- map_dfr(
  c("in_season", "out_season"),
  function(mode) {
    map_dfr(diagnosis_flags_fig, function(flag) {
      fit_single_outcome_panel_supp(build_diagnosis_events(flag), flag,
                               formula = ag_formula_3month,
                               season_mode = mode)
    }) %>%
      mutate(season_mode = mode)
  }
) %>%
  left_join(diagnosis_interpretation, by = c("outcome" = "diagnosis"))
print(diagnosis_panel_results_supp %>%
        select(season_mode, outcome, dx_group, exposure, HR, CI_low, CI_high,
               VE_pct, n_infants, n_events, fit_status),
      n = Inf)

panel_order <- c("RSV season", "Negative control")
arm_order   <- c("Neither", "Maternal RSV vaccine", "Nirsevimab")

fig3_fit_supp <- diagnosis_panel_results_supp %>%
  filter(str_starts(fit_status, "fit")) %>%
  mutate(
    panel_label = dplyr::recode(
      season_mode,
      in_season  = "RSV season",
      out_season = "Negative control"
    ),
    panel_label = factor(panel_label, levels = panel_order),

    exposure_label = case_when(
      exposure == "exposure_arma_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure == "exposure_armb_nirsevimab_only"   ~ "Nirsevimab",
      exposure == "exposure_armd_neither"            ~ "Neither"
    ),
    exposure_label = factor(exposure_label, levels = arm_order),

    row_kind = if_else(
      exposure_label == "Neither",
      "neither",
      "arm"
    ),

    outcome_label = dplyr::recode(
      outcome,
      !!!setNames(
        names(DIAGNOSIS_LABELS),
        unname(DIAGNOSIS_LABELS)
      )
    ),

    event_annotation = format_events_denom(
      n_events,
      n_infants,
      person_years
    ),

    ve_annotation = if_else(
      row_kind == "neither",
      "Reference",
      sprintf("%.1f (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high)
    )
  ) %>%
  filter(!is.na(exposure_label))

stopifnot(
  "No diagnoses had a usable fit in either panel -- nothing to plot" = nrow(fig3_fit_supp) > 0
)

diagnosis_catalog_all <- fig3_fit_supp %>%
  filter(row_kind == "arm") %>%
  distinct(outcome_label) %>%
  pull(outcome_label)

diagnosis_display_order <- fig3_fit_supp %>%
  filter(row_kind == "arm", outcome_label %in% diagnosis_catalog_all) %>%
  distinct(outcome_label, dx_group) %>%
  arrange(dx_group, outcome_label) %>%
  pull(outcome_label)
print(diagnosis_display_order)

FIG3_FONT_SIZE <- 11
LABEL_X      <- 0
LABEL_X_MIN  <- -0.15
INDENT_DX    <- 0.15
INDENT_ARM   <- 0.35
EVENTS_X     <- 4.2
TABLE_X_MAX  <- 5.8
VE_X_MIN     <- -0.15
VE_X_MAX     <- 2.2

transparent_bg <- theme(
  panel.background = element_blank(),
  plot.background  = element_rect(fill = "transparent", colour = NA)
)

build_season_figure <- function(season_lbl) {

  dx_in_group <- diagnosis_display_order

  scaffold <- map_dfr(dx_in_group, function(dx) {
    bind_rows(
      tibble(outcome_label = dx, exposure_label = NA_character_,
             row_type = "dx_header"),
      tibble(outcome_label = dx, exposure_label = arm_order,
             row_type = "data")
    )
  }) %>%
    mutate(
      row_kind = case_when(
        row_type == "data" & exposure_label == "Neither" ~ "neither",
        row_type == "data"                                ~ "arm",
        TRUE                                               ~ NA_character_
      ),
      row_id = row_number(), y_pos = -row_id,
      exposure_label_join = factor(exposure_label, levels = arm_order)
    )

  n_rows <- nrow(scaffold)
  Y_LIMITS <- c(-(n_rows + 0.5), -0.5)

  scaffold_data <- scaffold %>%
    filter(row_type == "data") %>%
    left_join(
      fig3_fit_supp %>%
        filter(panel_label == season_lbl, outcome_label %in% dx_in_group) %>%
        mutate(exposure_label_join = factor(exposure_label, levels = arm_order)) %>%
        select(outcome_label, exposure_label_join,
               VE_pct, VE_CI_low, VE_CI_high, event_annotation, ve_annotation),
      by = c("outcome_label", "exposure_label_join")
    )

  scaffold_full <- scaffold %>%
    left_join(scaffold_data %>% select(row_id, event_annotation, ve_annotation,
                                       VE_pct, VE_CI_low, VE_CI_high),
              by = "row_id") %>%
    mutate(
      ve_annotation = if_else(row_kind == "neither" & is.na(ve_annotation),
                              "Reference", ve_annotation)
    )

  shading <- scaffold_full %>%
    distinct(row_id, y_pos) %>%
    mutate(band = ifelse(row_id %% 2 == 0, "shade", "none")) %>%
    filter(band == "shade")

  p_bg <- ggplot(shading, aes(y = y_pos)) +
    geom_rect(
      aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5, xmin = -Inf, xmax = Inf),
      fill = "grey88", alpha = 0.5
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
    theme_void(base_size = FIG3_FONT_SIZE) +
    theme(
      plot.margin     = margin(t = 5.5, r = 0, b = 5.5, l = 0),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  p_table <- ggplot(scaffold_full, aes(y = y_pos)) +
    geom_text(
      data = filter(scaffold_full, row_type == "dx_header"),
      aes(x = LABEL_X + INDENT_DX, label = outcome_label), hjust = 0, fontface = "bold",
      size = FIG3_FONT_SIZE / .pt, color = "black"
    ) +
    geom_text(
      data = filter(scaffold_full, row_type == "data"),
      aes(x = LABEL_X + INDENT_ARM, label = exposure_label), hjust = 0,
      size = FIG3_FONT_SIZE / .pt, color = "black", na.rm = TRUE
    ) +
    geom_text(
      data = filter(scaffold_full, row_type == "data"),
      aes(x = EVENTS_X, label = event_annotation), hjust = 0,
      size = FIG3_FONT_SIZE / .pt, color = "black", na.rm = TRUE
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(LABEL_X_MIN, TABLE_X_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
    theme_void(base_size = FIG3_FONT_SIZE) +
    theme(plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5)) +
    transparent_bg

  table_header <- ggplot() +
    annotate("text", x = LABEL_X, y = 0.75, label = "Diagnosis and exposure arm", hjust = 0,
             fontface = "bold", size = FIG3_FONT_SIZE / .pt) +
    annotate("text", x = EVENTS_X, y = 0.75, label = fig_denom_header(), hjust = 0,
             fontface = "bold", size = FIG3_FONT_SIZE / .pt) +
    scale_x_continuous(limits = c(LABEL_X_MIN, TABLE_X_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 5.5, r = 0, b = 0, l = 5.5)) +
    transparent_bg

  p_forest <- ggplot(filter(scaffold_full, row_kind == "arm"), aes(x = VE_pct, y = y_pos)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = VE_CI_low, xmax = VE_CI_high), height = 0,
                   linewidth = 0.8, color = "black", na.rm = TRUE) +
    geom_point(size = 3.0, shape = 15, color = "black", na.rm = TRUE) +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 10),
      limits = c(-500, 76),
      breaks = c(-500, -150, -50, 0, 50),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
    labs(x = "Reduction in antibiotic prescriptions (%)", y = NULL) +
    theme_minimal(base_size = FIG3_FONT_SIZE) +
    theme(
      panel.grid    = element_blank(),
      axis.title.x  = element_text(size = FIG3_FONT_SIZE, color = "black"),
      axis.text.x   = element_text(size = FIG3_FONT_SIZE, color = "black"),
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.line.x   = element_line(color = "black", linewidth = 0.6),
      axis.ticks.x  = element_line(color = "black", linewidth = 0.6),
      legend.position = "none",
      plot.margin   = margin(t = 5.5, r = 0, b = 5.5, l = 0)
    ) +
    transparent_bg

  forest_header <- plot_spacer() + transparent_bg

  p_ve <- ggplot(scaffold_full, aes(y = y_pos)) +
    geom_text(
      data = filter(scaffold_full, row_type == "data"),
      aes(x = 0, label = ve_annotation), hjust = 0,
      size = FIG3_FONT_SIZE / .pt, color = "black", na.rm = TRUE
    ) +
    scale_x_continuous(limits = c(VE_X_MIN, VE_X_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
    theme_void(base_size = FIG3_FONT_SIZE) +
    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0)) +
    transparent_bg

  ve_header <- ggplot() +
    annotate("text", x = 0, y = 0.75, label = "VE, % (95% CI)", hjust = 0,
             fontface = "bold", size = FIG3_FONT_SIZE / .pt) +
    scale_x_continuous(limits = c(VE_X_MIN, VE_X_MAX), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 0)) +
    transparent_bg

  p_season_title <- ggplot() +
    annotate("text", x = 0, y = 1, label = as.character(season_lbl), hjust = 0,
             fontface = "bold", size = FIG3_FONT_SIZE / .pt, color = "black") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 5.5, r = 0, b = 2, l = 5.5)) +
    transparent_bg

  design <- c(
    area(1, 1, 1, 3),
    area(2, 1, 2, 1),
    area(2, 2, 2, 2),
    area(2, 3, 2, 3),
    area(3, 1, 3, 3),
    area(3, 1, 3, 1),
    area(3, 2, 3, 2),
    area(3, 3, 3, 3)
  )

  fig <- p_season_title +
    table_header + forest_header + ve_header +
    p_bg + p_table + p_forest + p_ve +
    plot_layout(design = design, heights = c(0.6, 1, n_rows * 0.9), widths = c(2.1, 1, 0.85))

  list(plot = fig, n_rows = n_rows)
}

seasons_present <- intersect(panel_order, unique(as.character(fig3_fit_supp$panel_label)))

season_figures <- map(seasons_present, build_season_figure)
names(season_figures) <- seasons_present

n_groups  <- length(season_figures)
ncol_grid <- n_groups

fig3_supp <- wrap_plots(map(season_figures, "plot"), ncol = ncol_grid, nrow = 1) +
  plot_annotation(
    title = "Supplementary Figure. Diagnosis-Specific VE, RSV Season and Negative Control",
    subtitle = "Never-both cohort, stratified by 3-month calendar period",
    
    theme = theme(
      plot.title    = element_text(size = FIG3_FONT_SIZE + 1, color = "black"),
      plot.subtitle = element_text(size = FIG3_FONT_SIZE, color = "black"),
      
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

print(fig3_supp)

max_rows  <- max(map_dbl(season_figures, "n_rows"))
fig_width <- 6.5 * ncol_grid

ggsave("figureS_diagnosis_inseason_negcontrol_0826.pdf", fig3_supp,
       width = fig_width + 3, height = 0.22 * max_rows + 2.1,
       units = "in", device = "pdf", limitsize = FALSE)
