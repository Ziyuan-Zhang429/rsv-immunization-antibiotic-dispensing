# Dispensing stacked
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(forcats)
})

if (!exists('cp_primary', inherits = TRUE) ||
    !exists('analytic_events', inherits = TRUE) ||
    !exists('cohort_primary', inherits = TRUE) ||
    !exists('clean_and_restrict', inherits = TRUE)) {
  source('R/primary/cox_primary.R')
}

out_dir <- 'output/antibiotic_dispensing_by_age'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_intervals <- clean_and_restrict(cp_primary$all_cause,
                                      restrict_to_season = FALSE)
all_event_days <- model_intervals %>%
  filter(event == 1L, tstop >= 1, tstop <= FOLLOWUP_MAX) %>%
  transmute(infant_study_id, abx_days = as.integer(tstop),
            exposure_arm = as.character(exposure_arm),
            gestage_imp = gestage_imp) %>%
  distinct(infant_study_id, abx_days, .keep_all = TRUE)

class_records <- analytic_events %>%
  filter(analyzable_event_primary == TRUE,
         !is.na(abx_days), abx_days >= 1, abx_days <= FOLLOWUP_MAX) %>%
  semi_join(cohort_primary %>% distinct(infant_study_id),
            by = 'infant_study_id') %>%
  mutate(
    antibiotic_group = case_when(
      abx_class == 'Aminoglycoside Antibiotics' ~ 'Aminoglycoside',
      abx_class == 'Vancomycin Antibiotics And Derivatives' ~ 'Vancomycin',
      abx_class == 'Penicillin Antibiotics' ~ 'Penicillin',
      abx_class %in% c('Cephalosporin Antibiotics - 1st Generation',
                       'Cephalosporin Antibiotics - 3rd Generation',
                       'Cephalosporin Anti - 4th Gen') ~ 'Cephalosporin',
      abx_class == 'Macrolide Antibiotics' ~ 'Macrolide',
      TRUE ~ 'Other antibiotic'
    )
  ) %>%
  distinct(infant_study_id, abx_days, antibiotic_group)

events <- class_records %>%
  inner_join(all_event_days, by = c('infant_study_id', 'abx_days'))

arm_labels <- c(
  d_neither = 'Neither',
  a_maternal_vax_only = 'Maternal RSVpreF',
  b_nirsevimab_only = 'Nirsevimab'
)
class_levels <- c('Aminoglycoside', 'Vancomycin', 'Penicillin',
                  'Cephalosporin', 'Macrolide', 'Other antibiotic')

plot_data <- bind_rows(
  events %>% mutate(gestational_group = 'All infants'),
  events %>% mutate(gestational_group = if_else(gestage_imp < 259,
                        'Preterm (32 to 36 weeks)', 'Term (≥37 weeks)'))
) %>%
  mutate(
    gestational_group = factor(gestational_group,
       levels = c('All infants', 'Preterm (32 to 36 weeks)', 'Term (≥37 weeks)')),
    exposure = factor(arm_labels[exposure_arm], levels = unname(arm_labels)),
    antibiotic_group = factor(antibiotic_group, levels = class_levels)
  )

# Plot the distribution of dispensing outcomes
make_stacked_plot <- function(max_day = 90L) {
  ggplot(filter(plot_data, abx_days <= max_day),
         aes(x = abx_days, fill = antibiotic_group)) +
    annotate('rect', xmin = 0.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
             fill = 'grey75', alpha = 0.22) +
    geom_bar(width = 1, color = NA) +
    facet_grid(gestational_group ~ exposure, scales = 'fixed', drop = FALSE) +
    scale_x_continuous(
      breaks = if (max_day <= 90L) c(1, 7, 14, 30, 60, 90)
               else c(1, 30, 60, 90, 120, 150, 180),
      limits = c(0.5, max_day + 0.5),
      expand = c(0, 0),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    coord_transform(y = scales::pseudo_log_trans(sigma = 1)) +
    scale_fill_brewer(palette = 'Set2', drop = FALSE) +
    labs(
      x = 'Day of life',
      y = 'Antibiotic-dispense days (count; pseudo-log scale)',
      fill = 'Antibiotic class'
    ) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(face = 'bold', size = 13),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 14, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          panel.spacing.x = grid::unit(1.5, 'lines'),
          panel.spacing.y = grid::unit(1.0, 'lines'),
          panel.grid.minor = element_blank(),
          legend.position = 'bottom') +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

p90 <- make_stacked_plot(90L)
p180 <- make_stacked_plot(180L)
ggsave(file.path(out_dir, 'antibiotic_dispensing_day1_90.pdf'),
       p90, width = 16, height = 11, units = 'in')
ggsave(file.path(out_dir, 'antibiotic_dispensing_day1_90.png'),
       p90, width = 16, height = 11, units = 'in', dpi = 300)
ggsave(file.path(out_dir, 'antibiotic_dispensing_day1_180.pdf'),
       p180, width = 18, height = 11, units = 'in')
ggsave(file.path(out_dir, 'antibiotic_dispensing_day1_180.png'),
       p180, width = 18, height = 11, units = 'in', dpi = 300)

period_labels <- c('Days 1–7', 'Days 8–28', 'Days 29–90', 'Days 91–180')
events_by_period <- plot_data %>%
  mutate(age_period = factor(case_when(
    abx_days <= 7 ~ 'Days 1–7',
    abx_days <= 28 ~ 'Days 8–28',
    abx_days <= 90 ~ 'Days 29–90',
    TRUE ~ 'Days 91–180'), levels = period_labels)) %>%
  count(gestational_group, exposure, age_period, antibiotic_group,
        name = 'n_dispense_days')
write_csv(events_by_period,
          file.path(out_dir, 'antibiotic_counts_by_period.csv'))

class_timing <- plot_data %>%
  filter(antibiotic_group %in% c('Aminoglycoside', 'Vancomycin')) %>%
  mutate(age_period = factor(case_when(
    abx_days <= 7 ~ 'Days 1–7',
    abx_days <= 28 ~ 'Days 8–28',
    abx_days <= 90 ~ 'Days 29–90',
    TRUE ~ 'Days 91–180'), levels = period_labels)) %>%
  count(gestational_group, exposure, antibiotic_group, age_period,
        name = 'n_dispense_days')
write_csv(class_timing,
          file.path(out_dir, 'aminoglycoside_vancomycin_timing.csv'))
