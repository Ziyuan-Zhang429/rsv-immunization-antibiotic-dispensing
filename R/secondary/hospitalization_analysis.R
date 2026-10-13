

suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr);library(readr);library(purrr);library(broom);library(survival)})
HOSP_SAS <- "kpsc_infant_hosp_10062026.sas7bdat"
LINK_FILE <- "official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv"
OUT_DIR <- "official_hospitalization_outputs"; dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
wcsv <- function(x,nm){write_csv(x,file.path(OUT_DIR,nm),na="");invisible(x)}
chr_id <- function(x) trimws(as.character(x))
for(x in c("cohort_primary","events_primary","build_cp","clean_and_restrict","fit_ag","ag_formula")) if(!exists(x,inherits=TRUE)) stop("Missing canonical object/function: ",x)
for(f in c(HOSP_SAS,LINK_FILE)) if(!file.exists(f)) stop("Missing: ",f)

link <- read_csv(LINK_FILE,show_col_types=FALSE) %>%
  mutate(across(c(infant_study_id_0113,infant_study_id_0616),chr_id))

valid_link <- link %>% filter(has_new_id,has_hospital_source_row,hospital_history_valid)
if(!nrow(valid_link)) stop("No hospital-valid officially linked infants")

h <- read_sas(HOSP_SAS) %>% transmute(
  infant_study_id_0616=chr_id(infant_study_id),
  admit_days=as.numeric(admit_days),disch_days=as.numeric(disch_days)) %>%
  inner_join(valid_link %>% select(infant_study_id_0113,infant_study_id_0616),by="infant_study_id_0616") %>%
  filter(!is.na(admit_days),!is.na(disch_days),admit_days>=0,disch_days>=admit_days) %>%
  distinct(infant_study_id_0113,admit_days,disch_days)

merge_one <- function(z){
  z <- arrange(z,admit_days,disch_days); a<-numeric(); b<-numeric()
  for(i in seq_len(nrow(z))){
    if(!length(b) || z$admit_days[i] > tail(b,1)){a<-c(a,z$admit_days[i]);b<-c(b,z$disch_days[i])}
    else b[length(b)]<-max(b[length(b)],z$disch_days[i])
  }
  tibble(admit_days=a,disch_days=b,episode_number=seq_along(a))
}
episodes <- bind_rows(lapply(split(h,h$infant_study_id_0113),merge_one),.id="infant_study_id_0113") %>%
  group_by(infant_study_id_0113) %>% arrange(admit_days,disch_days,.by_group=TRUE) %>%
  mutate(birth_day0=episode_number==1L & first(admit_days)==0,
         birth_day0_or_day1=episode_number==1L & first(admit_days)%in%c(0,1)) %>% ungroup()
wcsv(episodes %>% count(birth_day0,birth_day0_or_day1,name="episodes"),"01_episode_definition_QC.csv")

ev <- events_primary$all_cause %>% transmute(
  infant_study_id_0113=chr_id(infant_study_id),abx_days=as.numeric(abx_days)) %>%
  filter(abx_days>0,abx_days<=FOLLOWUP_MAX) %>% distinct() %>%
  semi_join(valid_link,by="infant_study_id_0113")
epi_split <- split(episodes,episodes$infant_study_id_0113)
class_one <- function(id,d){
  z<-epi_split[[id]]; if(is.null(z)) return(c(in_hospital=FALSE,birth_day0=FALSE,birth_day0_or_day1=FALSE))
  hit<-which(z$admit_days<=d & d<=z$disch_days)
  c(in_hospital=length(hit)>0,
    birth_day0=length(hit)>0 && any(z$birth_day0[hit]),
    birth_day0_or_day1=length(hit)>0 && any(z$birth_day0_or_day1[hit]))
}
flags<-t(mapply(class_one,ev$infant_study_id_0113,ev$abx_days))
ev<-bind_cols(ev,as_tibble(flags)) %>% mutate(category=case_when(
  birth_day0~"birth_stay_day0_definition",
  birth_day0_or_day1~"birth_stay_added_by_day1_definition",
  in_hospital~"later_hospital_overlap",TRUE~"outside_hospital"))
wcsv(ev %>% count(category,name="event_days"),"02_event_classification_QC.csv")
wcsv(ev %>% summarise(event_days=n(),during=sum(in_hospital),outside=sum(!in_hospital),
  birth_day0=sum(birth_day0),birth_day0_or_day1=sum(birth_day0_or_day1),
  added_by_day1=sum(birth_day0_or_day1 & !birth_day0)),"03_event_classification_totals.csv")
stopifnot(sum(ev$in_hospital)+sum(!ev$in_hospital)==nrow(ev))

sets<-list(
  H0_linked_baseline_all_cause=ev,
  H1_during_hospital=filter(ev,in_hospital),
  H2_outside_hospital=filter(ev,!in_hospital),
  H3_exclude_birth_stay_day0=filter(ev,!birth_day0),
  H4_exclude_birth_stay_day0_or_day1=filter(ev,!birth_day0_or_day1)
)
periods<-c("full_year","RSV_season","non_RSV_season")
valid_jan_ids<-valid_link$infant_study_id_0113
crude_formula<-Surv(tstart,tstop,event)~exposure_arm+strata(calendar_month)

fit_one_set<-function(e,label,period){
  e2<-e %>% transmute(infant_study_id=infant_study_id_0113,abx_days)

  if(is.numeric(cohort_primary$infant_study_id)) e2$infant_study_id<-as.numeric(e2$infant_study_id)
  cp<-build_cp(e2)
  md<-clean_and_restrict(cp,restrict_to_season=(period=="RSV_season"),
                         keep_out_of_season_only=(period=="non_RSV_season")) %>%
    filter(chr_id(infant_study_id)%in%valid_jan_ids)
  if(!nrow(md)) stop("No model data: ",label," / ",period)
  arms<-md %>% group_by(exposure_arm) %>% summarise(
    infants_contributing_time=n_distinct(infant_study_id),
    infants_with_events=n_distinct(infant_study_id[event==1]),event_days=sum(event),
    person_days=sum(tstop-tstart),.groups="drop") %>% mutate(analysis=label,period=period,.before=1)
  model_fit<-function(formula,model){
    fit<-tryCatch(fit_ag(md,formula=formula),error=function(e)e)
    if(inherits(fit,"error")) return(tibble(analysis=label,period,model,term=NA_character_,HR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,VE_pct=NA_real_,VE_CI_low=NA_real_,VE_CI_high=NA_real_,status=paste("ERROR",conditionMessage(fit))))
    tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(grepl("^exposure_arm",term)) %>% transmute(
      analysis=label,period,model,term,HR=estimate,CI_low=conf.low,CI_high=conf.high,
      VE_pct=100*(1-estimate),VE_CI_low=100*(1-conf.high),VE_CI_high=100*(1-conf.low),status="OK")
  }
  list(arms=arms,ve=bind_rows(model_fit(crude_formula,"calendar_stratified_unadjusted"),model_fit(ag_formula,"fully_adjusted")))
}
res<-list()
for(p in periods) for(nm in names(sets)) res[[paste(nm,p,sep="__")]]<-fit_one_set(sets[[nm]],nm,p)
arms<-bind_rows(map(res,"arms")); ve<-bind_rows(map(res,"ve"))
wcsv(arms,"04_events_infants_person_time_by_arm.csv");wcsv(ve,"05_VE_unadjusted_and_adjusted.csv")

term_map<-c(exposure_arma_maternal_vax_only="a_maternal_vax_only",exposure_armb_nirsevimab_only="b_nirsevimab_only")
main<-ve %>% filter(model=="fully_adjusted",term%in%names(term_map)) %>% mutate(exposure_arm=unname(term_map[term])) %>%
  left_join(arms,by=c("analysis","period","exposure_arm"))
ref<-arms %>% filter(as.character(exposure_arm)=="d_neither") %>% transmute(analysis,period,
  reference_infants=infants_contributing_time,reference_events=event_days,reference_person_days=person_days)
main<-main %>% left_join(ref,by=c("analysis","period")) %>% transmute(
  analysis,period,exposure_arm,reference_infants,reference_events,reference_person_days,
  exposed_infants=infants_contributing_time,exposed_events=event_days,exposed_person_days=person_days,
  VE_pct,VE_CI_low,VE_CI_high,status)
wcsv(main,"06_manuscript_style_adjusted_estimates.csv")
