

suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr);library(readr);library(stringr);library(purrr)})
DX_SAS <- "kpsc_infant_abx_dx_10062026.sas7bdat"
LINK_FILE <- "official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv"
OUT_DIR <- "official_D6_outputs";dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
wcsv<-function(x,nm){write_csv(x,file.path(OUT_DIR,nm),na="");invisible(x)}
chr_id<-function(x)trimws(as.character(x))
for(x in c("cohort_primary","analytic_events","FOLLOWUP_MAX")) if(!exists(x,inherits=TRUE)) stop("Missing canonical January object: ",x)
for(f in c(DX_SAS,LINK_FILE)) if(!file.exists(f)) stop("Missing: ",f)

classes<-c("Aminoglycoside Antibiotics","Vancomycin Antibiotics And Derivatives")
primary_ids<-cohort_primary %>% distinct(infant_study_id) %>% mutate(infant_study_id=chr_id(infant_study_id))
events<-analytic_events %>% filter(analyzable_event_primary==TRUE,!is.na(abx_days),abx_days>0,abx_days<=FOLLOWUP_MAX,abx_class%in%classes) %>%
  transmute(infant_study_id_0113=chr_id(infant_study_id),abx_class,abx_days=as.numeric(abx_days),across(any_of(c("asthma_allergy_yn","bronchitis_yn","fever_yn","influenza_yn","bacterial_infection_yn","other_infection_yn","sup_om_yn","nsup_om","pharyngitis_yn","pneumonia_yn","sinusitis_yn","UTI_yn","skin_cond_yn","URI_yn")))) %>%
  semi_join(primary_ids,by=c("infant_study_id_0113"="infant_study_id")) %>% distinct(infant_study_id_0113,abx_class,abx_days,.keep_all=TRUE) %>%
  mutate(age_group=cut(abx_days,c(0,7,30,90,180),labels=c("1-7","8-30","31-90","91-180"),right=TRUE))

age<-events %>% count(abx_class,age_group,name="n_event_days",.drop=FALSE) %>% group_by(abx_class) %>% mutate(total_event_days=sum(n_event_days),pct=100*n_event_days/total_event_days) %>% ungroup()
wcsv(age,"01_age_distribution_all_january_events.csv")

link<-read_csv(LINK_FILE,show_col_types=FALSE) %>% mutate(across(c(infant_study_id_0113,infant_study_id_0616),chr_id))
dx<-read_sas(DX_SAS) %>% transmute(infant_study_id_0616=chr_id(infant_study_id),dx_day=as.numeric(abx_dx_days),CODE=toupper(trimws(as.character(CODE)))) %>%
  filter(!is.na(dx_day),!is.na(CODE),CODE!="") %>% distinct()
dx_jan<-dx %>% inner_join(link %>% filter(has_new_id) %>% select(infant_study_id_0113,infant_study_id_0616),by="infant_study_id_0616")

make_links<-function(lo,hi,label){
  z<-events %>% inner_join(dx_jan,by="infant_study_id_0113",relationship="many-to-many") %>% mutate(relative_day=dx_day-abx_days) %>% filter(between(relative_day,lo,hi)) %>%
    distinct(infant_study_id_0113,abx_class,abx_days,age_group,dx_day,CODE,relative_day)
  keys<-z %>% distinct(infant_study_id_0113,abx_class,abx_days) %>% mutate(has_dx=TRUE)
  cov<-events %>% left_join(link %>% select(infant_study_id_0113,has_new_id,has_diagnosis_source_row),by="infant_study_id_0113") %>%
    left_join(keys,by=c("infant_study_id_0113","abx_class","abx_days")) %>% mutate(across(c(has_new_id,has_diagnosis_source_row,has_dx),~coalesce(.x,FALSE))) %>%
    group_by(abx_class,age_group) %>% summarise(window=label,n_all_event_days=n(),n_with_official_new_id=sum(has_new_id),n_with_dx_source_row=sum(has_diagnosis_source_row),n_with_dx_in_window=sum(has_dx),pct_with_dx_in_window=100*mean(has_dx),.groups="drop")
  list(links=z,coverage=cov)
}
main<-make_links(-4,0,"minus4_to_0"); pm7<-make_links(-7,7,"plusminus7")
wcsv(bind_rows(main$coverage,pm7$coverage),"02_diagnosis_coverage_by_age_and_window.csv")

cat_code<-function(code) case_when(
  grepl("^Z38",code)~"Birth-status (non-indication)",
  grepl("^Z05\\.?1",code)~"Observation for suspected neonatal infection (Z05.1)",
  grepl("^(P36|P37|P39|A4[0-1]|A49|B95|B96)",code)~"Infection-related diagnostic codes",
  grepl("^(J0[0-6]|J1[0-8]|J2[0-2]|P22|P23|P28)",code)~"Respiratory diagnostic codes",
  grepl("^P59",code)~"Neonatal jaundice",TRUE~"Other / not classified")
summary_window<-function(z,label){
  z<-z %>% mutate(dx_category=cat_code(CODE))
  den<-z %>% distinct(abx_class,infant_study_id_0113,abx_days) %>% count(abx_class,name="n_event_days_with_any_dx")
  z %>% distinct(abx_class,infant_study_id_0113,abx_days,dx_category) %>% count(abx_class,dx_category,name="n_event_days") %>%
    left_join(den,by="abx_class") %>% mutate(window=label,pct_of_events_with_any_dx=100*n_event_days/n_event_days_with_any_dx,.before=1)
}
wcsv(bind_rows(summary_window(main$links,"minus4_to_0"),summary_window(pm7$links,"plusminus7")),"03_broad_diagnosis_categories_by_window.csv")

top_codes<-function(z,label) z %>% filter(!grepl("^Z38",CODE)) %>% distinct(abx_class,infant_study_id_0113,abx_days,CODE) %>% count(abx_class,CODE,name="n_event_days") %>% group_by(abx_class) %>% slice_max(n_event_days,n=25,with_ties=FALSE) %>% ungroup() %>% mutate(window=label,.before=1)
wcsv(bind_rows(top_codes(main$links,"minus4_to_0"),top_codes(pm7$links,"plusminus7")),"04_top_nonbirth_ICD10_codes_by_window.csv")

first30_summary<-function(z,label){
  z %>% filter(abx_days<=30,!grepl("^Z38",CODE)) %>% distinct(abx_class,infant_study_id_0113,abx_days,CODE) %>% count(abx_class,CODE,name="n_first30_event_days") %>%
    group_by(abx_class) %>% arrange(desc(n_first30_event_days),.by_group=TRUE) %>% mutate(rank=row_number()) %>% filter(rank<=25) %>% ungroup() %>% mutate(window=label,.before=1)
}
wcsv(bind_rows(first30_summary(main$links,"minus4_to_0"),first30_summary(pm7$links,"plusminus7")),"05_first30_top_nonbirth_codes.csv")

first30_all_codes <- function(z, label) {
  z %>% filter(abx_days >= 1, abx_days <= 30) %>%
    distinct(abx_class, infant_study_id_0113, abx_days, CODE) %>%
    count(abx_class, CODE, name = "n_event_days") %>%
    mutate(window = label, .before = 1)
}
all_first30 <- bind_rows(
  first30_all_codes(main$links, "minus4_to_0"),
  first30_all_codes(pm7$links, "plusminus7")
)
wcsv(all_first30, "08_first30_all_ICD10_codes_long.csv")

all_denominators <- events %>% filter(abx_days >= 1, abx_days <= 30) %>%
  count(abx_class, name = "denominator_event_days")
all_coverage <- bind_rows(main$coverage, pm7$coverage) %>%
  filter(age_group %in% c("1-7", "8-30")) %>%
  group_by(abx_class, window) %>%
  summarise(n_with_dx_in_window = sum(n_with_dx_in_window), .groups = "drop")

for (abx in classes) {
  denom <- all_denominators %>% filter(abx_class == abx) %>% pull(denominator_event_days)
  code_table <- all_first30 %>% filter(abx_class == abx) %>%
    select(CODE, window, n_event_days) %>%
    pivot_wider(names_from = window, values_from = n_event_days, values_fill = 0) %>%
    mutate(
      `ICD10_CM_description` = NA_character_,
      denominator_event_days = denom,
      minus4_to_0_pct = round(100 * minus4_to_0 / denom, 1),
      plusminus7_pct = round(100 * plusminus7 / denom, 1)
    ) %>%
    arrange(desc(minus4_to_0), desc(plusminus7), CODE) %>%
    select(CODE, ICD10_CM_description, denominator_event_days,
           minus4_to_0, minus4_to_0_pct, plusminus7, plusminus7_pct)
  cov <- all_coverage %>% filter(abx_class == abx)
  n_main <- cov %>% filter(window == "minus4_to_0") %>% pull(n_with_dx_in_window)
  n_pm7 <- cov %>% filter(window == "plusminus7") %>% pull(n_with_dx_in_window)
  any_row <- tibble(
    CODE = "ANY_RECORDED_DIAGNOSIS", ICD10_CM_description = "Any recorded diagnosis",
    denominator_event_days = denom,
    minus4_to_0 = n_main, minus4_to_0_pct = round(100*n_main/denom,1),
    plusminus7 = n_pm7, plusminus7_pct = round(100*n_pm7/denom,1)
  )
  out <- bind_rows(code_table, any_row)
  name <- if (abx == classes[1]) "09_first30_aminoglycoside_ALL_codes.csv" else "10_first30_vancomycin_ALL_codes.csv"
  wcsv(out, name)
}

dx_flags<-intersect(c("asthma_allergy_yn","bronchitis_yn","fever_yn","influenza_yn","bacterial_infection_yn","other_infection_yn","sup_om_yn","nsup_om","pharyngitis_yn","pneumonia_yn","sinusitis_yn","UTI_yn","skin_cond_yn","URI_yn"),names(events))
jan14<-events %>% filter(abx_days<=30) %>% group_by(abx_class) %>% summarise(denominator_event_days=n(),across(all_of(dx_flags),~sum(.x==1,na.rm=TRUE)),.groups="drop") %>%
  pivot_longer(all_of(dx_flags),names_to="diagnosis_flag",values_to="n_event_days") %>% mutate(pct=100*n_event_days/denominator_event_days) %>% arrange(abx_class,desc(n_event_days))
wcsv(jan14,"06_first30_january_prespecified14_categories.csv")
none14<-events %>% filter(abx_days<=30) %>% mutate(any14=if_any(all_of(dx_flags),~coalesce(.x==1,FALSE))) %>% group_by(abx_class) %>% summarise(event_days=n(),none_of_14=sum(!any14),pct_none_of_14=100*mean(!any14),.groups="drop")
wcsv(none14,"07_first30_none_of_14_prespecified_categories.csv")
print(age);print(main$coverage);print(none14)
