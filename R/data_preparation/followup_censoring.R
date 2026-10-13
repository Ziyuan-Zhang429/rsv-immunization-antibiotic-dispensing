

suppressPackageStartupMessages({library(dplyr);library(readr);library(haven);library(purrr);library(survival)})
OUT_CENS <- 'censoring_revision_QC'; dir.create(OUT_CENS,showWarnings=FALSE,recursive=TRUE)
LINK_CENS <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
INF_CENS <- 'kpsc_infant_covariates_06162026.sas7bdat'
stopifnot(exists('cohort_primary'),exists('cohort_full'),exists('make_analysis_blocks'),exists('build_cp'),exists('events_primary'),exists('fit_ag'))
for(f in c(LINK_CENS,INF_CENS)) if(!file.exists(f)) stop('Missing: ',f)
id_chr <- function(x) trimws(as.character(x))
link <- read_csv(LINK_CENS,show_col_types=FALSE,col_types=cols(infant_study_id_0113=col_character(),infant_study_id_0616=col_character(),.default=col_guess())) %>%
  mutate(across(c(infant_study_id_0113,infant_study_id_0616),id_chr))
stopifnot(!anyDuplicated(link$infant_study_id_0113),!anyDuplicated(na.omit(link$infant_study_id_0616)))
new_cov <- read_sas(INF_CENS)
stopifnot(all(c('infant_study_id','disenroll_days') %in% names(new_cov)))
new_cov <- new_cov %>% transmute(infant_study_id_0616=id_chr(infant_study_id),disenroll_days=as.numeric(disenroll_days))
stopifnot(!anyDuplicated(new_cov$infant_study_id_0616))

master_cens <- cohort_primary %>% mutate(infant_study_id_0113=id_chr(infant_study_id)) %>%
  left_join(link %>% select(infant_study_id_0113,infant_study_id_0616),by='infant_study_id_0113') %>%
  left_join(new_cov,by='infant_study_id_0616') %>%
  mutate(has_new_id=!is.na(infant_study_id_0616),has_dis_days=!is.na(disenroll_days),
         followup_censor_days=if_else(has_dis_days,pmin(180,pmax(0,disenroll_days)),NA_real_))
stopifnot(nrow(master_cens)==nrow(cohort_primary),!anyDuplicated(master_cens$infant_study_id_0113))
saveRDS(master_cens,file.path(OUT_CENS,'january_primary_all_with_new_followup.rds'))
qc <- master_cens %>% summarise(original_n=n(),mapped=sum(has_new_id),with_dis_days=sum(has_dis_days),
  unlinked=sum(!has_new_id),linked_missing_dis_days=sum(has_new_id & !has_dis_days),
  disenroll_before_30=sum(disenroll_days<30,na.rm=TRUE),
  disenroll_before_180=sum(disenroll_days<180,na.rm=TRUE),
  implausibly_large=sum(disenroll_days>10000,na.rm=TRUE))
print(qc);write_csv(qc,file.path(OUT_CENS,'01_linkage_followup_qc.csv'))

CONFIRMED_DISENROLL_DEFINITION <- TRUE
if (!CONFIRMED_DISENROLL_DEFINITION) {
} else {
  linked_ids <- master_cens %>% filter(has_dis_days) %>% pull(infant_study_id_0113)
  censor_key <- master_cens %>% filter(has_dis_days) %>% select(infant_study_id_0113,followup_censor_days)
  attach_censor <- function(cohort) {
    cohort %>% mutate(infant_study_id_0113=id_chr(infant_study_id)) %>%
      inner_join(censor_key,by='infant_study_id_0113') %>%
      filter(followup_censor_days>0) %>% select(-infant_study_id_0113)
  }

  linked_cohort <- attach_censor(cohort_primary)
  original_linked_blocks <- make_analysis_blocks(linked_cohort,label='linked_uncensored')
  corrected_blocks <- original_linked_blocks
  corrected_blocks$base_surv <- linked_cohort %>%
    transmute(infant_study_id,tstart=0,tstop=followup_censor_days)
  stopifnot(nrow(corrected_blocks$base_surv)==nrow(linked_cohort),
            all(corrected_blocks$base_surv$tstop>0),
            all(corrected_blocks$base_surv$tstop<=180))
  run_variant <- function(blocks, label) {
    cp <- purrr::map(events_primary,~build_cp(.x,blocks=blocks))
    models <- list(in_season=purrr::map(cp,clean_and_restrict),
                   full_year=purrr::map(cp,~clean_and_restrict(.x,restrict_to_season=FALSE)),
                   out_of_season=purrr::map(cp,~clean_and_restrict(.x,keep_out_of_season_only=TRUE)))
    counts <- purrr::imap_dfr(models,function(lst,period) purrr::imap_dfr(lst,function(md,outcome) {
      tibble(variant=label,period,outcome,infants=n_distinct(md$infant_study_id),
             person_days=sum(md$tstop-md$tstart),events=sum(md$event))
    }))
    estimates <- purrr::imap_dfr(models,function(lst,period) purrr::imap_dfr(lst,function(md,outcome) {
      tryCatch({
        f <- fit_ag(md)
        tidy_ve(f,outcome) %>% mutate(variant=label,period,.before=1)
      },error=function(e) tibble(variant=label,period,outcome,exposure=NA_character_,fit_error=conditionMessage(e)))
    }))
    list(counts=counts,estimates=estimates,models=models)
  }
  unc <- run_variant(original_linked_blocks,'linked_original_180d')
  cor <- run_variant(corrected_blocks,'linked_censored')
  write_csv(bind_rows(unc$counts,cor$counts),file.path(OUT_CENS,'02_model_counts_comparison.csv'))
  write_csv(bind_rows(unc$estimates,cor$estimates),file.path(OUT_CENS,'03_ag_estimates_comparison.csv'))
  saveRDS(list(uncensored=unc,censored=cor),file.path(OUT_CENS,'linked_AG_comparison.rds'))
}
