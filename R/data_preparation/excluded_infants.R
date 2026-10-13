

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(haven)})
OUT <- 'revised_missingness_characteristics';dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
stopifnot(exists('analytic_infants'), 'has_complete_covariates' %in% names(analytic_infants))
id_clean <- function(x) {z<-trimws(as.character(x));z[z %in% c('','NA','NaN','<NA>')]<-NA_character_;z}
XWALK_SAS <- 'kpsc_study_id_xwalk_10092026.sas7bdat'
LINK_CSV <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
if(file.exists(XWALK_SAS)) {
  x<-read_sas(XWALK_SAS) %>% transmute(jan_id=id_clean(infant_study_id_0113),june_id=id_clean(infant_study_id_0616))
} else if(file.exists(LINK_CSV)) {
  x<-read_csv(LINK_CSV,col_types=cols(.default=col_character())) %>% transmute(jan_id=id_clean(infant_study_id_0113),june_id=id_clean(infant_study_id_0616))
} else stop('Official crosswalk SAS or QC CSV not found')
x<-x %>% filter(!is.na(jan_id))
stopifnot(!anyDuplicated(x$jan_id))
cohort<-analytic_infants %>% mutate(jan_id=id_clean(infant_study_id)) %>% left_join(x,by='jan_id') %>%
  mutate(has_linkage=!is.na(june_id), original_complete=has_complete_covariates %in% TRUE,
         revised_status=if_else(has_linkage & original_complete,'Included','Excluded'),
         exclusion_reason=case_when(original_complete & !has_linkage~'Unlinked only',
           !original_complete & has_linkage~'Original covariates incomplete only',
           !original_complete & !has_linkage~'Both missing covariates and linkage',
           TRUE~'Included'))
stopifnot(nrow(cohort)==47868,!anyDuplicated(cohort$jan_id))
print(cohort %>% count(revised_status,exclusion_reason))
write_csv(cohort %>% count(revised_status,exclusion_reason),file.path(OUT,'00_cohort_flow_counts.csv'))
stopifnot(sum(cohort$revised_status=='Included')==47219,sum(cohort$revised_status=='Excluded')==649)

field_labels<-c(exposure_group='RSV immunization group',gestage_imp='Gestational age at birth',
 apgar_1min='1-minute Apgar score',GENDER='Infant sex',race_cat_infant='Infant race/ethnicity',
 serious_newborn_comorbidity='Serious newborn comorbidity',age_index='Maternal age at delivery',
 race_cat_mom='Maternal race/ethnicity',bmi_prior='Maternal BMI (pre-pregnancy)',
 smoking='Maternal cigarette smoking',median_income_cat='Neighborhood median household income',
 charlson_Wt='Maternal weighted Charlson index',gbs_prior='Group B Streptococcus screening',
 flu_prior='Maternal influenza vaccination',tdap_prior='Maternal Tdap vaccination',
 OUTPT='Maternal outpatient visits, prior year',INPT='Maternal inpatient admissions, prior year',
 ED='Maternal emergency visits, prior year',tot_covid='Maternal COVID-19 vaccine doses',
 lb_prior='Prior live birth',rsv_valid='Maternal RSVpreF exposure',mab_history='Ever received nirsevimab')
vars<-intersect(names(field_labels),names(cohort))
if(length(setdiff(names(field_labels),vars))) invisible(NULL)
continuous<-intersect(c('gestage_imp','apgar_1min','age_index','bmi_prior','charlson_Wt','OUTPT','INPT','ED','tot_covid'),vars)
for(v in vars) if(inherits(cohort[[v]],'haven_labelled')) {
 if(v %in% continuous) cohort[[v]]<-as.numeric(cohort[[v]]) else cohort[[v]]<-haven::as_factor(cohort[[v]])
}
fmt_n<-function(n,d) sprintf('%d (%.1f%%)',n,100*n/d)
rows<-list();ranks<-list()
for(v in vars) {
 z<-cohort[[v]];lbl<-unname(field_labels[[v]]);a<-cohort$revised_status=='Included';b<-!a
 if(v %in% continuous && is.numeric(z)) {
  x<-as.numeric(z[a]);y<-as.numeric(z[b]);f<-function(t) if(all(is.na(t))) 'NA' else sprintf('%.2f (%.2f)',mean(t,na.rm=TRUE),sd(t,na.rm=TRUE))
  pool<-sqrt((var(x,na.rm=TRUE)+var(y,na.rm=TRUE))/2)
  smd<-if(is.finite(pool)&&pool>0) abs(mean(x,na.rm=TRUE)-mean(y,na.rm=TRUE))/pool else NA_real_
  rows[[length(rows)+1]]<-tibble(variable=v,characteristic=lbl,level='Mean (SD)',included=f(x),excluded=f(y),absolute_smd=smd)
  rows[[length(rows)+1]]<-tibble(variable=v,characteristic=lbl,level='Missing, n (%)',included=fmt_n(sum(is.na(x)),length(x)),excluded=fmt_n(sum(is.na(y)),length(y)),absolute_smd=NA_real_)
  ranks[[length(ranks)+1]]<-tibble(variable=v,characteristic=lbl,level='Continuous',absolute_smd=smd)
 } else {
  z<-as.character(z);z[is.na(z)|trimws(z)=='']<-'(Missing)'
  for(lev in sort(unique(z))) {
   p1<-mean(z[a]==lev);p0<-mean(z[b]==lev);den<-sqrt((p1*(1-p1)+p0*(1-p0))/2)
   smd<-if(den>0) abs(p1-p0)/den else 0
   rows[[length(rows)+1]]<-tibble(variable=v,characteristic=lbl,level=lev,included=fmt_n(sum(z[a]==lev),sum(a)),excluded=fmt_n(sum(z[b]==lev),sum(b)),absolute_smd=smd)
   ranks[[length(ranks)+1]]<-tibble(variable=v,characteristic=lbl,level=lev,absolute_smd=smd)
  }
 }
}
tab<-bind_rows(rows);detail<-bind_rows(ranks)
rank<-detail %>% group_by(variable,characteristic) %>% summarise(max_category_absolute_smd=if(all(is.na(absolute_smd))) NA_real_ else max(absolute_smd,na.rm=TRUE),.groups='drop') %>% arrange(desc(max_category_absolute_smd))
write_csv(tab,file.path(OUT,'01_included_vs_excluded_characteristics.csv'))
write_csv(detail,file.path(OUT,'02_category_level_SMDs.csv'))
write_csv(rank,file.path(OUT,'03_ranked_imbalances.csv'))

write_csv(cohort %>% summarise(eligible=n(),original_covariates_incomplete=sum(!original_complete),
 unavailable_linkage=sum(!has_linkage),both=sum(!original_complete & !has_linkage),
 total_excluded=sum(revised_status=='Excluded'),final_included=sum(revised_status=='Included')),
 file.path(OUT,'04_missingness_summary.csv'))
