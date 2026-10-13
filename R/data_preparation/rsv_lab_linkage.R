

suppressPackageStartupMessages({library(dplyr);library(readr);library(haven)})
OUT <- 'unified_jan_jun_oct'
LAB <- 'kpsc_infant_lb_outcomes_06162026.sas7bdat'
LINK <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
if (!file.exists(LAB) || !file.exists(LINK)) stop('Missing June lab SAS or official linkage QC file')
chr_id <- function(x) trimws(as.character(x))
link <- read_csv(LINK, show_col_types=FALSE, col_types=cols(infant_study_id_0113=col_character(),infant_study_id_0616=col_character(),.default=col_guess())) %>%
  filter(!is.na(infant_study_id_0616), infant_study_id_0616 != '') %>%
  select(infant_study_id_0113,infant_study_id_0616)
stopifnot(!anyDuplicated(link$infant_study_id_0113), !anyDuplicated(link$infant_study_id_0616))
lab <- read_sas(LAB)
stopifnot(all(c('infant_study_id','outcome_rsvlab_days','lab_type','result') %in% names(lab)))
lab <- lab %>% mutate(infant_study_id_0616=chr_id(infant_study_id))
linked <- lab %>% inner_join(link,by='infant_study_id_0616',relationship='many-to-one') %>%
  mutate(infant_study_id=infant_study_id_0113)

saveRDS(linked,file.path(OUT,'june_RSV_laboratory_long_mapped_to_january.rds'))
write_csv(tibble(metric=c('june_raw_lab_rows','june_raw_lab_infants','official_linked_lab_rows','official_linked_lab_infants','june_lab_rows_without_january_link','missing_lab_relative_day','linked_positive_or_detected_rows','linked_negative_rows'),
 value=c(nrow(lab),n_distinct(lab$infant_study_id_0616),nrow(linked),n_distinct(linked$infant_study_id_0113),nrow(lab)-nrow(linked),sum(is.na(linked$outcome_rsvlab_days)),sum(linked$result %in% c('Positive','Detected'),na.rm=TRUE),sum(linked$result %in% c('Negative','Not Detected'),na.rm=TRUE))),file.path(OUT,'june_RSV_lab_official_linkage_QC.csv'))

test_days <- linked %>% filter(!is.na(outcome_rsvlab_days),outcome_rsvlab_days>0,outcome_rsvlab_days<=180) %>%
  group_by(infant_study_id,outcome_rsvlab_days) %>% summarise(rsv_pos=any(result %in% c('Positive','Detected')),n_lab_rows=n(),.groups='drop') %>%
  mutate(rsv_neg=!rsv_pos)
saveRDS(test_days,file.path(OUT,'june_RSV_test_days_official.rds'))
