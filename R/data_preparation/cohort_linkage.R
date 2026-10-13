

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(readr); library(purrr)
})

XWALK_SAS <- "kpsc_study_id_xwalk_10092026.sas7bdat"
HOSP_SAS  <- "kpsc_infant_hosp_10062026.sas7bdat"
DX_SAS    <- "kpsc_infant_abx_dx_10062026.sas7bdat"
OUT_DIR   <- "official_crosswalk_QC"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

req_obj <- c("cohort_primary","events_primary","cp_primary","clean_and_restrict",
             "analytic_events","FOLLOWUP_MAX")
for (x in req_obj) if (!exists(x, inherits=TRUE)) stop("Missing canonical January object: ", x)
for (f in c(XWALK_SAS,HOSP_SAS,DX_SAS)) if (!file.exists(f)) stop("Missing input file: ", f)

wcsv <- function(x,nm) { write_csv(x,file.path(OUT_DIR,nm),na=""); invisible(x) }
chr_id <- function(x) trimws(as.character(x))
assert_cols <- function(x,cols,nm) {
  miss <- setdiff(cols,names(x)); if(length(miss)) stop(nm," missing: ",paste(miss,collapse=", "))
}

xw <- read_sas(XWALK_SAS)
assert_cols(xw,c("mom_study_id_0113","infant_study_id_0113",
                 "mom_study_id_0616","infant_study_id_0616"),"crosswalk")
xw <- xw %>% transmute(
  mom_study_id_0113=chr_id(mom_study_id_0113),
  infant_study_id_0113=chr_id(infant_study_id_0113),
  mom_study_id_0616=na_if(chr_id(mom_study_id_0616),"NA"),
  infant_study_id_0616=na_if(chr_id(infant_study_id_0616),"NA")
) %>% mutate(across(everything(), ~na_if(.x,"")))

xw_qc <- tibble(
  metric=c("rows","unique_old_infants","duplicate_old_infant_rows",
           "nonmissing_new_infant_rows","unique_nonmissing_new_infants",
           "duplicate_nonmissing_new_infant_rows","missing_new_infant_ids",
           "missing_new_mother_ids","discordant_new_id_missingness"),
  value=c(nrow(xw),n_distinct(xw$infant_study_id_0113),
          sum(duplicated(xw$infant_study_id_0113)),
          sum(!is.na(xw$infant_study_id_0616)),
          n_distinct(xw$infant_study_id_0616,na.rm=TRUE),
          sum(duplicated(xw$infant_study_id_0616) & !is.na(xw$infant_study_id_0616)),
          sum(is.na(xw$infant_study_id_0616)),sum(is.na(xw$mom_study_id_0616)),
          sum(xor(is.na(xw$infant_study_id_0616),is.na(xw$mom_study_id_0616))))
)
wcsv(xw_qc,"01_crosswalk_structure_QC.csv")
if (anyDuplicated(xw$infant_study_id_0113)) stop("Official crosswalk old infant IDs are not unique")
if (anyDuplicated(xw$infant_study_id_0616[!is.na(xw$infant_study_id_0616)])) stop("Official crosswalk new infant IDs are not unique")

cohort_jan <- cohort_primary %>% mutate(
  infant_study_id_0113=chr_id(infant_study_id),
  mom_study_id_0113_cohort=chr_id(mom_study_id)
)
cohort_map <- cohort_jan %>% left_join(xw,by="infant_study_id_0113") %>% mutate(
  in_official_crosswalk=!is.na(mom_study_id_0113),
  has_new_id=!is.na(infant_study_id_0616),
  mother_id_0113_agrees=is.na(mom_study_id_0113) | mom_study_id_0113_cohort==mom_study_id_0113
)
if (any(!cohort_map$mother_id_0113_agrees)) stop("January mother ID disagrees with official crosswalk for mapped infant(s)")

md_full <- clean_and_restrict(cp_primary$all_cause, restrict_to_season=FALSE,
                              keep_out_of_season_only=FALSE)
cc_ids <- md_full %>% distinct(infant_study_id) %>% transmute(infant_study_id_0113=chr_id(infant_study_id))

coverage_by_population <- bind_rows(
  cohort_map %>% summarise(population="cohort_primary_pre_complete_case",
    n=n(),in_crosswalk=sum(in_official_crosswalk),with_new_id=sum(has_new_id),
    pct_with_new_id=100*mean(has_new_id)),
  cohort_map %>% semi_join(cc_ids,by="infant_study_id_0113") %>% summarise(
    population="canonical_full_year_complete_case",
    n=n(),in_crosswalk=sum(in_official_crosswalk),with_new_id=sum(has_new_id),
    pct_with_new_id=100*mean(has_new_id))
)
wcsv(coverage_by_population,"02_january_population_mapping_coverage.csv")

exposure_summary <- cohort_map %>% mutate(
  mapping_status=if_else(has_new_id,"mapped_new_id","no_new_id"),
  maternal_rsv=as.integer(rsv_valid %in% TRUE),
  ever_nirsevimab=as.integer(mab_history %in% TRUE)
) %>% group_by(mapping_status) %>% summarise(
  infants=n(),mothers=n_distinct(mom_study_id_0113_cohort),
  maternal_rsv_n=sum(maternal_rsv,na.rm=TRUE),
  maternal_rsv_pct=100*mean(maternal_rsv,na.rm=TRUE),
  ever_nirsevimab_n=sum(ever_nirsevimab,na.rm=TRUE),
  ever_nirsevimab_pct=100*mean(ever_nirsevimab,na.rm=TRUE),.groups="drop")
wcsv(exposure_summary,"03_mapped_vs_unmapped_exposure_QC.csv")

jan_events <- events_primary$all_cause %>% transmute(
  infant_study_id_0113=chr_id(infant_study_id),abx_days=as.numeric(abx_days)) %>%
  filter(abx_days>0,abx_days<=FOLLOWUP_MAX) %>% distinct()
event_map_qc <- cohort_map %>% select(infant_study_id_0113,has_new_id) %>%
  left_join(jan_events %>% count(infant_study_id_0113,name="event_days"),by="infant_study_id_0113") %>%
  mutate(event_days=coalesce(event_days,0L)) %>% group_by(has_new_id) %>% summarise(
    infants=n(),infants_with_event=sum(event_days>0),event_days=sum(event_days),
    event_days_per_100_infants=100*event_days/infants,.groups="drop")
wcsv(event_map_qc,"04_mapped_vs_unmapped_antibiotic_event_QC.csv")

hosp <- read_sas(HOSP_SAS); dx <- read_sas(DX_SAS)
assert_cols(hosp,c("infant_study_id","admit_days","disch_days"),"hospital file")
assert_cols(dx,c("infant_study_id","abx_dx_days","CODE"),"diagnosis file")
hosp <- hosp %>% mutate(infant_study_id_0616=chr_id(infant_study_id),
                        admit_days=as.numeric(admit_days),disch_days=as.numeric(disch_days))
dx <- dx %>% mutate(infant_study_id_0616=chr_id(infant_study_id),
                    abx_dx_days=as.numeric(abx_dx_days),CODE=toupper(trimws(as.character(CODE))))

mapped <- cohort_map %>% filter(has_new_id) %>% select(infant_study_id_0113,infant_study_id_0616)
hosp_ids <- hosp %>% distinct(infant_study_id_0616)
dx_ids <- dx %>% distinct(infant_study_id_0616)
source_cov <- mapped %>% mutate(
  has_hospital_source_row=infant_study_id_0616 %in% hosp_ids$infant_study_id_0616,
  has_diagnosis_source_row=infant_study_id_0616 %in% dx_ids$infant_study_id_0616
) %>% summarise(
  mapped_january_infants=n(),
  with_hospital_source_row=sum(has_hospital_source_row),
  pct_with_hospital_source_row=100*mean(has_hospital_source_row),
  with_diagnosis_source_row=sum(has_diagnosis_source_row),
  pct_with_diagnosis_source_row=100*mean(has_diagnosis_source_row)
)
wcsv(source_cov,"05_october_source_coverage.csv")

hosp_qc <- hosp %>% semi_join(mapped,by="infant_study_id_0616") %>% summarise(
  rows=n(),unique_infants=n_distinct(infant_study_id_0616),
  missing_admit=sum(is.na(admit_days)),missing_discharge=sum(is.na(disch_days)),
  negative_admit=sum(admit_days<0,na.rm=TRUE),
  discharge_before_admit=sum(disch_days<admit_days,na.rm=TRUE),
  exact_duplicate_intervals=n()-nrow(distinct(cur_data(),infant_study_id_0616,admit_days,disch_days)),
  day0_admission_rows=sum(admit_days==0,na.rm=TRUE),day1_admission_rows=sum(admit_days==1,na.rm=TRUE)
)
wcsv(hosp_qc,"06_hospital_interval_QC.csv")

bad_hosp <- hosp %>% semi_join(mapped,by="infant_study_id_0616") %>%
  filter(is.na(admit_days)|is.na(disch_days)|admit_days<0|disch_days<admit_days) %>%
  count(infant_study_id_0616,name="bad_interval_rows")
wcsv(bad_hosp,"07_hospital_infants_with_invalid_intervals.csv")

dx_qc <- dx %>% semi_join(mapped,by="infant_study_id_0616") %>% summarise(
  rows=n(),unique_infants=n_distinct(infant_study_id_0616),
  missing_dx_day=sum(is.na(abx_dx_days)),missing_or_blank_code=sum(is.na(CODE)|CODE==""),
  exact_duplicate_rows=n()-nrow(distinct(cur_data(),infant_study_id_0616,abx_dx_days,CODE)),
  z051_rows=sum(grepl("^Z05\\.?1",CODE),na.rm=TRUE)
)
wcsv(dx_qc,"08_diagnosis_extract_QC.csv")

bad_ids <- bad_hosp$infant_study_id_0616
linked_ids <- cohort_map %>% transmute(
  infant_study_id_0113,mom_study_id_0113=mom_study_id_0113_cohort,
  infant_study_id_0616,mom_study_id_0616,
  has_new_id,
  has_hospital_source_row=infant_study_id_0616 %in% hosp_ids$infant_study_id_0616,
  hospital_history_valid=has_hospital_source_row & !(infant_study_id_0616 %in% bad_ids),
  has_diagnosis_source_row=infant_study_id_0616 %in% dx_ids$infant_study_id_0616,
  in_canonical_full_year_complete_case=infant_study_id_0113 %in% cc_ids$infant_study_id_0113
)
wcsv(linked_ids,"09_official_linkage_status_by_january_infant.csv")
print(coverage_by_population)
print(source_cov)
print(hosp_qc)
print(dx_qc)
