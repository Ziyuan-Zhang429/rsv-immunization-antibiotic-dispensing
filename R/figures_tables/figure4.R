

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)

ADJUSTED_CSV <- "adjusted_absolute_effects.csv"
CRUDE_CSV    <- "crude_absolute_effects.csv"

for (f in c(ADJUSTED_CSV, CRUDE_CSV)) {
  if (!file.exists(f)) stop("Could not find ", f, ". Put it in the working directory or update the path at the top of this script.")
}

adjusted <- read.csv(ADJUSTED_CSV, stringsAsFactors = FALSE) %>%
  transmute(
    panel = period,
    exposure,
    followup = age_window,
    n_recipients_contributing = recipients_standardized,
    estimate,
    ci_low = ci_lower,
    ci_high = ci_upper
  )

crude <- read.csv(CRUDE_CSV, stringsAsFactors = FALSE) %>%
  transmute(
    panel,
    exposure,
    followup,
    crude_n_recipients = n_recipients_contributing,
    crude_estimate = estimate,
    crude_ci_low = ci_low,
    crude_ci_high = ci_high
  )

normalize_followup <- function(x) {
  x <- gsub("[–—]", "-", x)
  x <- gsub("^1-90 days$", "1-90 days", x)
  x <- gsub("^1-180 days$", "1-180 days", x)
  x
}
normalize_panel <- function(x) {
  x <- trimws(x)
  ifelse(x %in% c("RSV season", "RSV seasons"), "RSV season", x)
}

adjusted <- adjusted %>%
  mutate(panel = normalize_panel(panel), followup = normalize_followup(followup))
crude <- crude %>%
  mutate(panel = normalize_panel(panel), followup = normalize_followup(followup))

figure4_data <- adjusted %>%
  left_join(crude, by = c("panel", "exposure", "followup")) %>%
  mutate(
    panel = factor(panel, levels = c("Full year", "RSV season")),
    followup = factor(followup, levels = c("1-90 days", "1-180 days"))
  ) %>%
  arrange(panel, factor(exposure, levels = c("Maternal RSV vaccine", "Nirsevimab")), followup)

if (nrow(figure4_data) != 8L) stop("Expected 8 Figure 4 rows after joining the two result tables; found ", nrow(figure4_data), ".")
if (any(is.na(figure4_data$crude_estimate))) stop("Some adjusted rows did not match crude_absolute_effects.csv. Check panel/exposure/followup labels.")
if (any(figure4_data$n_recipients_contributing != figure4_data$crude_n_recipients)) {
  warning("Recipient counts differ between adjusted and crude files. Check that both analyses used the same analytic population.")
}
print(tibble::as_tibble(figure4_data), n = Inf)

make_courses_averted_plot <- function(df) {
  FONT <- 11; PAD <- 2
  panel_order <- c("Full year", "RSV season")
  exposure_order <- c("Maternal RSV vaccine", "Nirsevimab")
  headers <- crossing(panel = panel_order, exposure = exposure_order) %>%
    mutate(followup = NA_character_, n_recipients_contributing = NA_integer_, estimate = NA_real_,
           ci_low = NA_real_, ci_high = NA_real_, crude_n_recipients = NA_integer_,
           crude_estimate = NA_real_, crude_ci_low = NA_real_, crude_ci_high = NA_real_,
           row_type = "exposure")
  panel_headers <- tibble(panel = panel_order, exposure = NA_character_, followup = NA_character_,
                          n_recipients_contributing = NA_integer_, estimate = NA_real_, ci_low = NA_real_,
                          ci_high = NA_real_, crude_n_recipients = NA_integer_, crude_estimate = NA_real_,
                          crude_ci_low = NA_real_, crude_ci_high = NA_real_, row_type = "panel")
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
      adjusted_annotation = if_else(
        row_type == "data", sprintf("%.2f (%.2f, %.2f)", estimate, ci_low, ci_high), NA_character_),
      crude_annotation = if_else(
        row_type == "data", sprintf("%.2f (%.2f, %.2f)", crude_estimate, crude_ci_low, crude_ci_high), NA_character_)
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

  forest_adjusted <- plot_df %>%
    filter(row_type == "data") %>%
    transmute(y_pos, est = estimate, low = ci_low, high = ci_high)

  p_forest <- ggplot(forest_adjusted, aes(x = est, y = y_pos)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = .5) +
    geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = .7, color = "black") +
    geom_point(shape = 15, size = 2.8, color = "black") +
    scale_x_continuous(limits = c(-0.5, 5), breaks = pretty(c(0, 5), n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) +
    labs(x = "Antibiotic courses averted per 100 intervention recipients", y = NULL, shape = NULL) +
    theme_minimal(base_size = FONT) +
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = FONT, color = "black"),
          axis.text.x = element_text(size = FONT, color = "black"), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.x = element_line(color = "black", linewidth = .6),
          axis.ticks.x = element_line(color = "black", linewidth = .6), legend.position = "none",
          plot.margin = margin(5.5, 8, 5.5, 8)) + transparent

  forest_header <- plot_spacer() + transparent

  p_est <- ggplot(plot_df, aes(y = y_pos)) +
    geom_text(data = filter(plot_df, row_type == "data"),
              aes(x = 0, label = crude_annotation), hjust = 0, size = FONT / .pt) +
    geom_text(data = filter(plot_df, row_type == "data"),
              aes(x = 2.55, label = adjusted_annotation), hjust = 0, size = FONT / .pt) +
    scale_x_continuous(limits = c(0, 5.25), expand = c(0, 0)) +
    scale_y_continuous(limits = Y, expand = c(0, 0)) + theme_void(base_size = FONT) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, PAD)) + transparent

  est_header <- ggplot() +
    annotate("text", x = 2.625, y = 1.55,
             label = "Antibiotic courses averted per 100 recipients (95% CI)",
             hjust = .5, fontface = "bold", size = FONT / .pt) +
    annotate("text", x = 0, y = .65, label = "Unadjusted", hjust = 0,
             fontface = "bold", size = FONT / .pt) +
    annotate("text", x = 2.55, y = .65, label = "Adjusted", hjust = 0,
             fontface = "bold", size = FONT / .pt) +
    scale_x_continuous(limits = c(0, 5.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.3, 1.95), expand = c(0, 0)) + theme_void() +
    theme(plot.margin = margin(5.5, 5.5, 0, PAD)) + transparent

  design <- c(area(1,1,1,1), area(1,2,1,2), area(1,3,1,3), area(2,1,2,3),
              area(2,1,2,1), area(2,2,2,2), area(2,3,2,3))
  table_header + est_header + forest_header + p_bg + p_table + p_est + p_forest +
    plot_layout(design = design, heights = c(2.2, n_rows), widths = c(1.20, 1.55, .85)) +
    plot_annotation(
      title = "Antibiotic Courses Averted per 100 Intervention Recipients",
      
      theme = theme(plot.title = element_text(size = FONT),
                    
                    plot.background = element_rect(fill = "white", colour = NA)))
}

fig4 <- make_courses_averted_plot(figure4_data)
print(fig4)
ggsave("figure4_courses_averted_crude_adjusted.pdf", fig4,
       width = 10.2, height = 5.5, units = "in", device = "pdf")
