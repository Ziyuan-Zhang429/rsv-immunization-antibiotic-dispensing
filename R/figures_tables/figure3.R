

library(patchwork)

stopifnot(exists("arm_counts"), exists("check_arm_counts"),
          exists("format_events_denom"), exists("CANONICAL_N"))

validate_figure2_model_data <- function(md, label) {
  stopifnot(all(md$tstart >= 0), all(md$tstop <= FOLLOWUP_MAX),
            all(md$tstop > md$tstart))
  duplicated_events <- md %>% filter(event == 1L) %>%
    count(infant_study_id, tstop, name = "n") %>% filter(n > 1)
  if (nrow(duplicated_events)) {
    stop(label, ": multiple antibiotic events found for the same infant-day.")
  }
  invisible(TRUE)
}
validate_figure2_model_data(model_data_primary$all_cause, "Primary RSV-season model")
validate_figure2_model_data(model_data_fullyear$all_cause, "Full-year model")
validate_figure2_model_data(model_data_out_of_season$all_cause, "Negative-control model")

FIG2_PANELS <- list(
  a_rsv_season       = list(md = model_data_primary$all_cause,       label = "RSV seasons"),
  b_full_year        = list(md = model_data_fullyear$all_cause,      label = "Full year"),
  c_negative_control = list(md = model_data_out_of_season$all_cause, label = "Non-RSV seasons")
)

fig2_counts <- imap_dfr(FIG2_PANELS, function(p, panel) {
  arm_counts(p$md, p$label) %>%
    mutate(panel = panel, panel_label = p$label) %>%
    select(panel, panel_label, exposure, exposure_arm,
           n_infants_contributing, person_days, person_years, n_events)
})

fig2_panel_n <- imap_dfr(FIG2_PANELS, ~ tibble(
  panel = .y, panel_label = .x$label,
  n_panel = n_distinct(.x$md$infant_study_id)
))

for (nm in names(FIG2_PANELS)) {
  invisible(check_arm_counts(FIG2_PANELS[[nm]]$md, FIG2_PANELS[[nm]]$label))
}

ref_exposure_code <- "exposure_armd_neither"

CRUDE_RELATIVE_CSV <- "crude_relative_effects.csv"
if (!file.exists(CRUDE_RELATIVE_CSV)) {
  stop("Could not find ", CRUDE_RELATIVE_CSV,
       ". Put crude_relative_effects.csv in the working directory or update CRUDE_RELATIVE_CSV.")
}

crude_relative <- read.csv(CRUDE_RELATIVE_CSV, stringsAsFactors = FALSE) %>%
  filter(outcome == "all_cause") %>%
  mutate(
    panel = recode(period,
                   "RSV seasons" = "a_rsv_season",
                   "Full year" = "b_full_year",
                   "Non-RSV seasons" = "c_negative_control")
  ) %>%
  select(panel, exposure,
         crude_HR = HR, crude_CI_low = CI_low, crude_CI_high = CI_high,
         crude_VE_pct = VE_pct, crude_VE_CI_low = VE_CI_low,
         crude_VE_CI_high = VE_CI_high)

stopifnot(nrow(crude_relative) == 6L)

fig2_ve <- bind_rows(
  primary_ve_table       %>% filter(outcome == "all_cause") %>% mutate(panel = "a_rsv_season"),
  sens_ve_table_fullyear %>% filter(outcome == "all_cause") %>% mutate(panel = "b_full_year"),
  out_of_season_ve_table %>% filter(outcome == "all_cause") %>% mutate(panel = "c_negative_control")
) %>%
  filter(exposure %in% c("exposure_arma_maternal_vax_only", "exposure_armb_nirsevimab_only")) %>%
  select(panel, exposure, HR, CI_low, CI_high, VE_pct, VE_CI_low, VE_CI_high) %>%
  left_join(crude_relative, by = c("panel", "exposure")) %>%
  left_join(fig2_counts, by = c("panel", "exposure")) %>%
  mutate(
    exposure_label = case_when(
      exposure == "exposure_arma_maternal_vax_only" ~ "Maternal RSV vaccine",
      exposure == "exposure_armb_nirsevimab_only"   ~ "Nirsevimab"
    ),
    row_type = "data",
    event_annotation = format_events_denom(n_events, n_infants_contributing, person_years),
    adjusted_ve_annotation = sprintf("%.1f (%.1f, %.1f)", VE_pct, VE_CI_low, VE_CI_high),
    crude_ve_annotation = sprintf("%.1f (%.1f, %.1f)", crude_VE_pct, crude_VE_CI_low, crude_VE_CI_high),
    ve_annotation = adjusted_ve_annotation
  )

stopifnot(!any(is.na(fig2_ve$n_events)), !any(is.na(fig2_ve$person_years)))

fig2_neither <- fig2_counts %>%
  filter(exposure == ref_exposure_code) %>%
  mutate(
    exposure_label = "Neither",
    row_type = "neither",
    HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
    crude_VE_pct = NA_real_, crude_VE_CI_low = NA_real_, crude_VE_CI_high = NA_real_,
    adjusted_ve_annotation = NA_character_, crude_ve_annotation = NA_character_,
    event_annotation = format_events_denom(n_events, n_infants_contributing, person_years),
    ve_annotation = NA_character_
  )

fig2_headers <- fig2_panel_n %>%
  mutate(
    row_type = "header",
    exposure_label = sprintf("%s", panel_label, scales::comma(n_panel)),
    HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    VE_pct = NA_real_, VE_CI_low = NA_real_, VE_CI_high = NA_real_,
    crude_VE_pct = NA_real_, crude_VE_CI_low = NA_real_, crude_VE_CI_high = NA_real_,
    adjusted_ve_annotation = NA_character_, crude_ve_annotation = NA_character_,
    event_annotation = NA_character_, ve_annotation = NA_character_
  )

panel_order <- c("RSV seasons", "Full year", "Non-RSV seasons")
arm_order   <- c("Maternal RSV vaccine", "Nirsevimab")

fig2_data <- map_dfr(panel_order, function(pl) {
  bind_rows(
    fig2_headers %>% filter(panel_label == pl),
    fig2_neither %>% filter(panel_label == pl),
    fig2_ve %>% filter(panel_label == pl) %>%
      mutate(exposure_label = factor(exposure_label, levels = arm_order)) %>%
      arrange(exposure_label) %>%
      mutate(exposure_label = as.character(exposure_label))
  )
}) %>%
  mutate(row_id = row_number(), y_pos = -row_id)

n_rows <- nrow(fig2_data)

fig2_shading <- fig2_data %>%
  distinct(row_id, y_pos) %>%
  mutate(band = ifelse(row_id %% 2 == 0, "shade", "none")) %>%
  filter(band == "shade")
print(fig2_data %>%
        select(row_id, panel_label, row_type, exposure_label,
               n_infants_contributing, person_years, n_events,
               event_annotation, ve_annotation),
      n = Inf)
print(CANONICAL_N$complete_case %>% select(exposure, n_infants_cohort))

FIG2_FONT_SIZE <- 11
Y_LIMITS <- c(-(n_rows + 0.5), -0.5)

LABEL_X <- 0
LABEL_X_MIN <- -0.4
INDENT_ARM <- 0.15
EVENTS_X <- 2.65
TABLE_X_MAX <- 4
PAD <- 2

transparent_bg <- theme(
  panel.background = element_blank(),
  plot.background  = element_rect(fill = "transparent", colour = NA)
)

p_table <- ggplot(fig2_data, aes(y = y_pos)) +
  geom_text(data = filter(fig2_data, row_type == "header"),
            aes(x = LABEL_X, label = exposure_label), hjust = 0, fontface = "bold",
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  geom_text(data = filter(fig2_data, row_type %in% c("data", "neither")),
            aes(x = LABEL_X + INDENT_ARM, label = exposure_label), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  geom_text(data = filter(fig2_data, row_type %in% c("data", "neither")),
            aes(x = EVENTS_X, label = event_annotation), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  scale_x_continuous(limits = c(LABEL_X_MIN, TABLE_X_MAX), expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
  theme_void(base_size = FIG2_FONT_SIZE) +
  theme(plot.margin = margin(5.5, PAD, 5.5, 5.5)) +
  transparent_bg

table_header <- ggplot() +
  annotate("text", x = LABEL_X, y = 0.6, label = "Exposure arm",
           hjust = 0, fontface = "bold", size = FIG2_FONT_SIZE / .pt) +
  annotate("text", x = EVENTS_X, y = 0.6, label = fig_denom_header(),
           hjust = 0, fontface = "bold", size = FIG2_FONT_SIZE / .pt) +
  scale_x_continuous(limits = c(LABEL_X_MIN, TABLE_X_MAX), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5.5, PAD, 0, 5.5)) +
  transparent_bg

fig2_forest_adjusted <- fig2_data %>%
  filter(row_type == "data") %>%
  transmute(y_pos, est = VE_pct, low = VE_CI_low, high = VE_CI_high)

p_forest <- ggplot(fig2_forest_adjusted, aes(x = est, y = y_pos)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = low, xmax = high),
                 height = 0, linewidth = 0.7, color = "black") +
  geom_point(shape = 15, size = 2.8, color = "black") +
  scale_x_continuous(limits = c(-25, 50), breaks = seq(-20, 40, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
  labs(x = "Reduction in all-cause antibiotic prescriptions (%)", y = NULL, shape = NULL) +
  theme_minimal(base_size = FIG2_FONT_SIZE) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = FIG2_FONT_SIZE, color = "black"),
    axis.text.x = element_text(size = FIG2_FONT_SIZE, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.margin = margin(5.5, PAD, 5.5, PAD)
  ) +
  transparent_bg

forest_header <- plot_spacer() + transparent_bg

p_ve <- ggplot(fig2_data, aes(y = y_pos)) +
  geom_text(data = filter(fig2_data, row_type == "data"),
            aes(x = 0, label = crude_ve_annotation), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  geom_text(data = filter(fig2_data, row_type == "data"),
            aes(x = 2.25, label = adjusted_ve_annotation), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  geom_text(data = filter(fig2_data, row_type == "neither"),
            aes(x = 0, label = "Reference"), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  geom_text(data = filter(fig2_data, row_type == "neither"),
            aes(x = 2.25, label = "Reference"), hjust = 0,
            size = FIG2_FONT_SIZE / .pt, color = "black") +
  scale_x_continuous(limits = c(0, 4.65), expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
  theme_void(base_size = FIG2_FONT_SIZE) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, PAD)) +
  transparent_bg

ve_header <- ggplot() +
  annotate("text", x = 2.325, y = 1.55,
           label = "Treatment effect, % (95% CI)",
           hjust = 0.5, fontface = "bold", size = FIG2_FONT_SIZE / .pt) +
  annotate("text", x = 0, y = 0.65, label = "Unadjusted",
           hjust = 0, fontface = "bold", size = FIG2_FONT_SIZE / .pt) +
  annotate("text", x = 2.25, y = 0.65, label = "Adjusted",
           hjust = 0, fontface = "bold", size = FIG2_FONT_SIZE / .pt) +
  scale_x_continuous(limits = c(0, 4.65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 1.95), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5.5, 5.5, 0, PAD)) +
  transparent_bg

fig2_design <- c(
  area(1,1,1,1), area(1,2,1,2), area(1,3,1,3),
  area(2,1,2,3), area(2,1,2,1), area(2,2,2,2), area(2,3,2,3)
)

p_bg <- ggplot(fig2_shading, aes(y = y_pos)) +
  geom_rect(aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5, xmin = -Inf, xmax = Inf),
            fill = "grey88", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = Y_LIMITS, expand = c(0, 0)) +
  theme_void(base_size = FIG2_FONT_SIZE) +
  theme(plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 0),
        plot.background = element_rect(fill = "white", colour = NA))

fig2 <- table_header + ve_header + forest_header +
  p_bg + p_table + p_ve + p_forest +
  plot_layout(design = fig2_design, heights = c(2.2, n_rows),
              widths = c(1.25, 1.55, 0.95)) +
  plot_annotation(
    title = "Effectiveness of RSV Immunization Against Antibiotic Use",
    subtitle = "Forest plot (right) shows adjusted estimates; both models reported numerically",
    
    theme = theme(
      plot.title = element_text(size = FIG2_FONT_SIZE),
      plot.subtitle = element_text(size = FIG2_FONT_SIZE),
      
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

print(fig2)

ggsave("figure3_main_ve_forest_crude_adjusted.pdf", fig2, width = 10.5, height = 6.1,
       units = "in", device = "pdf")

library(gt)

stopifnot(
  exists("arm_counts"),
  exists("check_arm_counts"),
  exists("ARM_COUNT_NOTE")
)

TABLE2_PANELS <- list(
  rsv_season = list(
    md = model_data_primary$all_cause,
    label = "RSV season"
  ),
  full_year = list(
    md = model_data_fullyear$all_cause,
    label = "Full year"
  ),
  out_season = list(
    md = model_data_out_of_season$all_cause,
    label = "Out of RSV season"
  )
)

table2_counts <- imap_dfr(TABLE2_PANELS, function(p, period_code) {
  arm_counts(p$md, p$label) %>%
    mutate(
      period_code = period_code,
      period = p$label
    ) %>%
    select(
      period_code,
      period,
      exposure,
      n_infants = n_infants_contributing,
      person_days,
      n_events
    )
})

for (nm in names(TABLE2_PANELS)) {
  invisible(
    check_arm_counts(
      TABLE2_PANELS[[nm]]$md,
      TABLE2_PANELS[[nm]]$label
    )
  )
}

table2_ve <- bind_rows(
  primary_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period_code = "rsv_season"),

  sens_ve_table_fullyear %>%
    filter(outcome == "all_cause") %>%
    mutate(period_code = "full_year"),

  out_of_season_ve_table %>%
    filter(outcome == "all_cause") %>%
    mutate(period_code = "out_season")
) %>%
  filter(
    exposure %in% c(
      "exposure_arma_maternal_vax_only",
      "exposure_armb_nirsevimab_only"
    )
  ) %>%
  select(
    period_code,
    exposure,
    HR,
    CI_low,
    CI_high,
    VE_pct,
    VE_CI_low,
    VE_CI_high
  )

table2_neither <- table2_counts %>%
  filter(exposure == "exposure_armd_neither") %>%
  transmute(
    period_code,
    period,
    exposure,
    n_infants,
    person_days,
    n_events,
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    VE_pct = NA_real_,
    VE_CI_low = NA_real_,
    VE_CI_high = NA_real_
  )

table2_exposed <- table2_ve %>%
  left_join(
    table2_counts,
    by = c("period_code", "exposure")
  )

table2 <- bind_rows(
  table2_exposed,
  table2_neither
) %>%
  mutate(
    `Analysis period` = factor(
      period,
      levels = c(
        "RSV season",
        "Full year",
        "Out of RSV season"
      )
    ),

    Exposure = factor(
      case_when(
        exposure == "exposure_armd_neither" ~
          "Neither",

        exposure == "exposure_arma_maternal_vax_only" ~
          "Maternal RSV vaccine",

        exposure == "exposure_armb_nirsevimab_only" ~
          "Nirsevimab"
      ),
      levels = c(
        "Neither",
        "Maternal RSV vaccine",
        "Nirsevimab"
      )
    ),

    `HR (95% CI)` = if_else(
      Exposure == "Neither",
      "Reference",
      sprintf(
        "%.2f (%.2f, %.2f)",
        HR,
        CI_low,
        CI_high
      )
    ),

    `VE, % (95% CI)` = if_else(
      Exposure == "Neither",
      "Reference",
      sprintf(
        "%.1f (%.1f, %.1f)",
        VE_pct,
        VE_CI_low,
        VE_CI_high
      )
    )
  ) %>%
  arrange(
    `Analysis period`,
    Exposure
  ) %>%
  select(
    `Analysis period`,
    Exposure,
    `Infants contributing person-time` = n_infants,
    `Person-days` = person_days,
    `Antibiotic events` = n_events,
    `HR (95% CI)`,
    `VE, % (95% CI)`
  )
print(table2, n = Inf)

table2_gt <- table2 %>%
  gt(
    groupname_col = "Analysis period"
  ) %>%
  fmt_number(
    columns = c(
      `Infants contributing person-time`,
      `Person-days`,
      `Antibiotic events`
    ),
    decimals = 0,
    use_seps = TRUE
  ) %>%
  tab_header(
    title = md(
      "**Table 2. Effectiveness of maternal RSV vaccination and nirsevimab against all-cause antibiotic use during infancy**"
    )
  ) %>%
  tab_source_note(
    source_note = md(
      paste0(
        "VE = vaccine/monoclonal-antibody effectiveness; ",
        "HR = hazard ratio. Antibiotic events were defined as unique ",
        "infant-days with systemic antibiotic use during days 1–180. ",
        "Models were adjusted for prespecified covariates and stratified ",
        "by calendar month, with infant age as the Cox time scale. ",
        "Neither served as the reference exposure group. ",
        "The out-of-RSV-season analysis served as a negative-control analysis. ",
        ARM_COUNT_NOTE
      )
    )
  )

print(table2_gt)

gtsave(
  table2_gt,
  "table2_main_all_cause_effectiveness_sensitivity_window.docx"
)
