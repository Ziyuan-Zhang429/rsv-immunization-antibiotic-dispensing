

library(dplyr)
library(tibble)
stopifnot(exists('events_primary'), exists('analytic_events'),
          exists('cohort_primary'), exists('model_data_fullyear'),
          exists('apply_washout'), exists('FOLLOWUP_MAX'))
keys <- c('infant_study_id','abx_days')
normalize <- function(x) x %>%
  filter(!is.na(abx_days), abx_days > 0, abx_days <= FOLLOWUP_MAX) %>%
  distinct(infant_study_id, abx_days)
all_events <- normalize(events_primary$all_cause)
primary_ids <- cohort_primary %>% distinct(infant_study_id)
model_ids <- model_data_fullyear$all_cause %>% distinct(infant_study_id)
primary_events <- all_events %>% semi_join(primary_ids, by='infant_study_id')
model_eligible_events <- primary_events %>% semi_join(model_ids, by='infant_study_id')

rebuilt_primary <- analytic_events %>%
  filter(analyzable_event_primary == TRUE) %>%
  normalize() %>% semi_join(primary_ids, by='infant_study_id')
count_row <- function(stage, x) tibble(stage=stage, event_days=nrow(x),
  infants_with_event=n_distinct(x$infant_study_id))
stages <- bind_rows(
  count_row('Eligible cohort, all event days', all_events),
  count_row('Primary cohort (exclude ever-both)', primary_events),
  count_row('Primary + full-year model-eligible IDs', model_eligible_events),
  count_row('Independent primary rebuild from analytic_events', rebuilt_primary))
stages$removed_since_previous <- c(NA_integer_, -diff(stages$event_days[1:3]), NA_integer_)
print(stages, n=Inf)
if (exists('model_data_primary')) invisible(NULL)
if (exists('n_primary_all_cause_events')) invisible(NULL)
if (exists('n_primary_all_cause_events_eligible')) invisible(NULL)

for (nm in c('all_events','primary_events','model_eligible_events')) {
  x <- get(nm)
}
write.csv(stages,'antibiotic_event_count_audit.csv',row.names=FALSE)

