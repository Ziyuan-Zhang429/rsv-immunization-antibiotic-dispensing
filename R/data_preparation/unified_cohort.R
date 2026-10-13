

suppressPackageStartupMessages({library(dplyr);library(readr);library(haven);library(tidyr)})
OUT <- 'unified_jan_jun_oct'; dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
chr_id <- function(x) trimws(as.character(x))
required <- c('cohort_primary','events_primary','cp_primary','clean_and_restrict')
for (nm in required) if (!exists(nm,inherits=TRUE)) stop('Source canonical January script first: ',nm)
LINK <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
HOSP <- 'kpsc_infant_hosp_10062026.sas7bdat'
DX <- 'kpsc_infant_abx_dx_10062026.sas7bdat'
for (f in c(LINK,HOSP,DX)) if (!file.exists(f)) stop('Missing ',f)
link <- read_csv(LINK,show_col_types=FALSE,col_types=cols(infant_study_id_0113=col_character(),mom_study_id_0113=col_character(),infant_study_id_0616=col_character(),mom_study_id_0616=col_character(),.default=col_guess()))
stopifnot(!anyDuplicated(link$infant_study_id_0113),!anyDuplicated(na.omit(link$infant_study_id_0616)))
cohort <- cohort_primary %>% mutate(infant_study_id_0113=chr_id(infant_study_id),mom_study_id_0113_check=chr_id(mom_study_id)) %>%
  left_join(link,by='infant_study_id_0113')
stopifnot(nrow(cohort)==nrow(cohort_primary),all(cohort$mom_study_id_0113_check==cohort$mom_study_id_0113))
cohort <- cohort %>% select(-mom_study_id_0113_check)

write_rds(cohort,file.path(OUT,'cohort_spine_january_canonical.rds'))
write_csv(cohort %>% select(infant_study_id_0113,mom_study_id_0113,infant_study_id_0616,mom_study_id_0616,has_new_id,has_hospital_source_row,hospital_history_valid,has_diagnosis_source_row,in_canonical_full_year_complete_case),file.path(OUT,'cohort_linkage_flags.csv'))
map <- link %>% filter(has_new_id) %>% select(infant_study_id_0113,infant_study_id_0616)
hosp <- read_sas(HOSP) %>% mutate(infant_study_id_0616=chr_id(infant_study_id)) %>%
  inner_join(map,by='infant_study_id_0616')
dx <- read_sas(DX) %>% mutate(infant_study_id_0616=chr_id(infant_study_id)) %>%
  inner_join(map,by='infant_study_id_0616')

write_rds(hosp,file.path(OUT,'october_hospital_long_mapped_to_january.rds'))
write_rds(dx,file.path(OUT,'october_diagnosis_long_mapped_to_january.rds'))
write_rds(events_primary,file.path(OUT,'january_primary_antibiotic_event_sets.rds'))

attach_june <- function(path,table_name,infant_col='infant_study_id') {
  z <- read_sas(path) %>% mutate(infant_study_id_0616=chr_id(.data[[infant_col]])) %>%
    inner_join(map,by='infant_study_id_0616')
  saveRDS(z,file.path(OUT,paste0('june_',table_name,'_long_mapped_to_january.rds')))
  invisible(z)
}

qc <- tibble(metric=c('january_cohort_rows','mapped_new_id','hospital_source_covered','diagnosis_source_covered','hospital_long_rows','diagnosis_long_rows','hospital_unique_january_infants','diagnosis_unique_january_infants'),
 value=c(nrow(cohort),sum(cohort$has_new_id),sum(cohort$has_hospital_source_row),sum(cohort$has_diagnosis_source_row),nrow(hosp),nrow(dx),n_distinct(hosp$infant_study_id_0113),n_distinct(dx$infant_study_id_0113)))
write_csv(qc,file.path(OUT,'unified_merge_QC.csv'))
print(qc)
