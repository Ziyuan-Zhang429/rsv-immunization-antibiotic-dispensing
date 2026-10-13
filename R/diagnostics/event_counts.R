

library(dplyr)
library(tidyr)

stopifnot(exists('model_data_primary'), exists('model_data_fullyear'),
          exists('model_data_out_of_season'))

summarize_abx_arm <- function(md, analysis_name) {
  required <- c('infant_study_id', 'exposure_arm', 'tstart', 'tstop', 'event')
  stopifnot(all(required %in% names(md)))
  stopifnot(all(md$tstop > md$tstart), !anyNA(md$event),
            all(md$event %in% c(0, 1)))

  out <- md %>%
    group_by(exposure_arm) %>%
    summarise(
      n_infants_contributing = n_distinct(infant_study_id),
      person_days = sum(tstop - tstart),
      antibiotic_events = sum(event),
      infants_with_event_in_arm = n_distinct(infant_study_id[event == 1]),
      .groups = 'drop'
    ) %>%
    mutate(
      analysis = analysis_name,
      exposure = recode(as.character(exposure_arm),
        'd_neither' = 'Neither (unexposed person-time)',
        'a_maternal_vax_only' = 'Maternal RSVpreF',
        'b_nirsevimab_only' = 'Nirsevimab',
        'c_both' = 'Both products'),
      person_years = person_days / 365.25,
      events_per_100_person_years = if_else(
        person_years > 0, 100 * antibiotic_events / person_years, NA_real_
      )
    ) %>%
    select(analysis, exposure, n_infants_contributing,
           person_days, person_years, antibiotic_events,
           infants_with_event_in_arm, events_per_100_person_years)

  stopifnot(sum(out$antibiotic_events) == sum(md$event),
            isTRUE(all.equal(sum(out$person_days), sum(md$tstop - md$tstart))))
  out
}

reviewer5_counts <- bind_rows(
  summarize_abx_arm(model_data_primary$all_cause, 'RSV season (primary)'),
  summarize_abx_arm(model_data_fullyear$all_cause, 'Full year'),
  summarize_abx_arm(model_data_out_of_season$all_cause, 'Outside RSV season')
)

print(reviewer5_counts, n = Inf, width = Inf)
reviewer5_totals <- reviewer5_counts %>%
  group_by(analysis) %>%
  summarise(person_days = sum(person_days),
            antibiotic_events = sum(antibiotic_events), .groups = 'drop')
print(reviewer5_totals, n = Inf)
write.csv(reviewer5_counts, 'reviewer5_comment3_counts.csv', row.names = FALSE)
write.csv(reviewer5_totals, 'reviewer5_comment3_totals.csv', row.names = FALSE)
