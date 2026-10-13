

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})
OUT <- 'official_RSV_neonatal_testing_outputs'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
wcsv <- function(x, name) { write_csv(x, file.path(OUT,name), na=''); invisible(x) }
chr_id <- function(x) trimws(as.character(x))
required <- c('cohort_primary','analytic_events','FOLLOWUP_MAX')
for (nm in required) if (!exists(nm,inherits=TRUE)) stop('Missing canonical January object: ',nm)
LAB_RDS <- 'unified_jan_jun_oct/june_RSV_laboratory_long_mapped_to_january.rds'
LINK_CSV <- 'official_crosswalk_QC/09_official_linkage_status_by_january_infant.csv'
for (f in c(LAB_RDS, LINK_CSV)) if (!file.exists(f)) stop('Missing prerequisite: ',f)

SEASONS <- data.frame(
  start=as.Date(c('2023-10-01','2024-11-01')),
  end=as.Date(c('2024-02-29','2025-03-31'))
)
is_rsv_season <- function(d) {
  vapply(as.Date(d), function(x) {
    if (is.na(x)) return(NA)
    any(x >= SEASONS$start & x <= SEASONS$end)
  }, logical(1))
}

cohort <- cohort_primary %>% mutate(infant_id=chr_id(infant_study_id))
stopifnot(!anyDuplicated(cohort$infant_id))
DATE_CANDIDATES <- c('index_date','birth_date','dob','date_of_birth','infant_dob','birth_dt')
found <- intersect(DATE_CANDIDATES, names(cohort))
if (!length(found)) stop('Cannot identify infant birth/index date in cohort_primary. Available columns: ',paste(names(cohort),collapse=', '), '\nSet BIRTH_DATE_COLUMN explicitly in this script after confirming its meaning.')
BIRTH_DATE_COLUMN <- found[1]
cohort <- cohort %>% mutate(birth_date=as.Date(.data[[BIRTH_DATE_COLUMN]]))
if (anyNA(cohort$birth_date)) stop('Birth date contains NA/unparseable values; confirm BIRTH_DATE_COLUMN and source format.')

expected_start <- as.Date("2023-06-01")
expected_end   <- as.Date("2025-06-30")

date_margin <- 7L

outside_margin <- cohort$birth_date < (expected_start - date_margin) |
  cohort$birth_date > (expected_end + date_margin)

if (any(outside_margin)) {
  stop("Some index dates fall more than 7 days outside the expected cohort range.")
}
cohort <- cohort %>% select(infant_id,birth_date)

link <- read_csv(LINK_CSV, show_col_types=FALSE,
  col_types=cols(infant_study_id_0113=col_character(),infant_study_id_0616=col_character(),.default=col_guess())) %>%
  transmute(infant_id=chr_id(infant_study_id_0113), newer_id=chr_id(infant_study_id_0616)) %>%
  mutate(newer_id=na_if(newer_id,'NA'),newer_id=na_if(newer_id,''))
stopifnot(!anyDuplicated(link$infant_id), !anyDuplicated(na.omit(link$newer_id)))
cohort <- cohort %>% left_join(link,by='infant_id') %>% mutate(has_official_new_id=!is.na(newer_id))

CLASSES <- c('Aminoglycoside Antibiotics','Vancomycin Antibiotics And Derivatives')
events <- analytic_events %>%
  filter(analyzable_event_primary==TRUE, !is.na(abx_days), abx_days>=1,abx_days<=30,
         abx_class %in% CLASSES) %>%
  transmute(infant_id=chr_id(infant_study_id),abx_class,abx_day=as.numeric(abx_days)) %>%
  distinct() %>% inner_join(cohort,by='infant_id') %>%
  mutate(abx_date=birth_date+abx_day,
         season=if_else(is_rsv_season(abx_date),'RSV season','Outside RSV season'))
if (!nrow(events)) stop('No January first-30-day aminoglycoside/vancomycin events.')

lab_raw <- readRDS(LAB_RDS)
needed <- c('infant_study_id_0113','outcome_rsvlab_days','result')
if (!all(needed %in% names(lab_raw))) stop('Mapped lab RDS missing: ',paste(setdiff(needed,names(lab_raw)),collapse=', '))
lab <- lab_raw %>% transmute(
  infant_id=chr_id(infant_study_id_0113),
  test_day=as.numeric(outcome_rsvlab_days),
  result_text=tolower(trimws(as.character(result)))) %>%
  mutate(status=case_when(result_text %in% c('positive','detected') ~ 'positive',
                          result_text %in% c('negative','not detected') ~ 'negative',
                          TRUE ~ 'other'))
if (any(lab$status=='other')) {
  wcsv(lab %>% filter(status=='other') %>% count(result_text,sort=TRUE),'00_unclassified_lab_results.csv')
  warning('Unclassified RSV lab results found; see 00_unclassified_lab_results.csv. These do not enter positivity denominators.')
}
lab <- lab %>% filter(!is.na(test_day),test_day>=0,test_day<=180) %>%
  inner_join(cohort %>% select(infant_id,birth_date),by='infant_id') %>%
  mutate(test_date=birth_date+test_day,
         season=if_else(is_rsv_season(test_date),'RSV season','Outside RSV season'),
         week_of_life=floor(test_day/7)+1L)

test_days <- lab %>% group_by(infant_id,test_day,test_date,season,week_of_life) %>%
  summarise(any_positive=any(status=='positive'),any_negative=any(status=='negative'),
            has_interpretable=any(status %in% c('positive','negative')),
            n_raw_rows=n(),.groups='drop') %>%
  mutate(test_result=case_when(any_positive~'positive',any_negative~'negative',TRUE~'unclassified'))
wcsv(test_days %>% count(test_result,name='test_days'),'00_test_day_result_QC.csv')

recipient_groups <- events %>% group_by(abx_class,season,infant_id) %>%
  summarise(n_abx_event_days=n(),first_abx_day=min(abx_day),.groups='drop') %>%
  left_join(cohort %>% select(infant_id,has_official_new_id),by='infant_id')
infant_lab <- test_days %>% filter(test_day>=1,test_day<=30,has_interpretable) %>%
  group_by(infant_id) %>% summarise(any_test_first30=TRUE,
      any_positive_first30=any(any_positive),.groups='drop')
recipient_groups <- recipient_groups %>% left_join(infant_lab,by='infant_id') %>%
  mutate(any_test_first30=coalesce(any_test_first30,FALSE),
         any_positive_first30=coalesce(any_positive_first30,FALSE))
A_infant <- recipient_groups %>% group_by(abx_class,season) %>% summarise(
  n_recipients=n(),n_mapped=sum(has_official_new_id),n_unmapped=sum(!has_official_new_id),
  n_mapped_tested_first30=sum(any_test_first30 & has_official_new_id),
  n_mapped_positive_first30=sum(any_positive_first30 & has_official_new_id),
  .groups='drop') %>% mutate(
  pct_mapped_tested_first30=if_else(n_mapped>0,100*n_mapped_tested_first30/n_mapped,NA_real_),
  n_mapped_without_recorded_test=n_mapped-n_mapped_tested_first30,
  pct_mapped_without_recorded_test=if_else(n_mapped>0,100*n_mapped_without_recorded_test/n_mapped,NA_real_),
  positivity_among_tested_pct=if_else(n_mapped_tested_first30>0,100*n_mapped_positive_first30/n_mapped_tested_first30,NA_real_))
wcsv(A_infant,'01A_antibiotic_recipients_testing_by_season.csv')

matched <- events %>% select(infant_id,abx_class,abx_day,season,has_official_new_id) %>%
  left_join(test_days %>% filter(has_interpretable) %>%
              select(infant_id,test_day,any_positive,test_result),
            by='infant_id',relationship='many-to-many') %>%
  mutate(delta=test_day-abx_day)
A_event <- matched %>% group_by(infant_id,abx_class,abx_day,season,has_official_new_id) %>%
  summarise(any_pm7=any(!is.na(delta)&between(delta,-7,7)),
    any_before=any(!is.na(delta)&between(delta,-7,0)),
    any_after=any(!is.na(delta)&between(delta,0,7)),
    any_positive_pm7=any(!is.na(delta)&between(delta,-7,7)&any_positive,na.rm=TRUE),
    any_positive_before=any(!is.na(delta)&between(delta,-7,0)&any_positive,na.rm=TRUE),
    .groups='drop')
A_event_summary <- A_event %>% group_by(abx_class,season) %>% summarise(
  n_event_days=n(),n_mapped_event_days=sum(has_official_new_id),
  n_mapped_tested_pm7=sum(any_pm7 & has_official_new_id),
  n_mapped_tested_before_or_same=sum(any_before & has_official_new_id),
  n_mapped_tested_after_or_same=sum(any_after & has_official_new_id),
  n_mapped_positive_pm7=sum(any_positive_pm7 & has_official_new_id),
  n_mapped_positive_before_or_same=sum(any_positive_before & has_official_new_id),
  .groups='drop') %>% mutate(across(starts_with('n_mapped_tested_'),
    ~if_else(n_mapped_event_days>0,100*.x/n_mapped_event_days,NA_real_),.names='pct_{.col}'))
wcsv(A_event_summary,'01B_antibiotic_event_testing_windows_by_season.csv')

MAX_WEEK <- 8L
weekly_lab <- test_days %>% filter(week_of_life<=MAX_WEEK,has_interpretable) %>%
  group_by(infant_id,week_of_life,season) %>% summarise(
    any_positive=any(any_positive),n_test_days=n(),.groups='drop')

week_grid <- tidyr::crossing(cohort %>% filter(has_official_new_id) %>% select(infant_id,birth_date),
                             week_of_life=seq_len(MAX_WEEK)) %>%
  mutate(week_start=birth_date+(week_of_life-1L)*7,
         week_end=week_start+6L)

season_overlap <- function(a,b) any(a<=SEASONS$end & b>=SEASONS$start)
week_grid <- week_grid %>% mutate(has_season=mapply(season_overlap,week_start,week_end))
den <- bind_rows(
  week_grid %>% filter(has_season) %>% count(week_of_life,name='n_cohort_infants') %>% mutate(season='RSV season'),
  week_grid %>% filter(!has_season | (week_start<min(SEASONS$start) | week_end>max(SEASONS$end))) %>%

    count(week_of_life,name='n_cohort_infants') %>% mutate(season='Outside RSV season'))

outside_any <- function(a,b) {
  days <- seq(a,b,by='day'); any(!is_rsv_season(days))
}

outside_flag <- mapply(outside_any,week_grid$week_start,week_grid$week_end)
den_out <- week_grid[outside_flag,] %>% count(week_of_life,name='n_cohort_infants') %>% mutate(season='Outside RSV season')
den <- bind_rows(week_grid %>% filter(has_season) %>% count(week_of_life,name='n_cohort_infants') %>% mutate(season='RSV season'),den_out)
weekly <- weekly_lab %>% group_by(week_of_life,season) %>% summarise(
  n_tested_infants=n(),n_positive_infants=sum(any_positive),.groups='drop') %>%
  right_join(den,by=c('week_of_life','season')) %>%
  mutate(n_tested_infants=coalesce(n_tested_infants,0L),
         n_positive_infants=coalesce(n_positive_infants,0L),
         n_tested_negative=n_tested_infants-n_positive_infants,
         positivity_pct=if_else(n_tested_infants>0,100*n_positive_infants/n_tested_infants,NA_real_),
         pct_cohort_tested=if_else(n_cohort_infants>0,100*n_tested_infants/n_cohort_infants,NA_real_)) %>%
  arrange(season,week_of_life)
wcsv(weekly,'02_weekly_RSV_testing_positivity_full_cohort.csv')

DX_RDS <- 'unified_jan_jun_oct/october_diagnosis_long_mapped_to_january.rds'
if (file.exists(DX_RDS)) {
  dx0 <- readRDS(DX_RDS)
  dx_required <- c('infant_study_id_0113','abx_dx_days','CODE')
  if (all(dx_required %in% names(dx0))) {
    z05 <- dx0 %>% transmute(infant_id=chr_id(infant_study_id_0113),
        dx_day=as.numeric(abx_dx_days),code=toupper(trimws(as.character(CODE)))) %>%
      filter(code=='Z05.1',!is.na(dx_day)) %>% distinct()

    selected_events <- events %>%
      filter(abx_class %in% CLASSES[1:2])

    z05_match <- selected_events %>%
      inner_join(z05, by = "infant_id",
                 relationship = "many-to-many") %>%
      filter(between(dx_day - abx_day, -4, 0)) %>%
      distinct(infant_id, abx_day) %>%
      mutate(z05_flag = TRUE)

    z05_events <- selected_events %>%
      left_join(z05_match, by = c("infant_id", "abx_day")) %>%
      mutate(z05_flag = coalesce(z05_flag, FALSE)) %>%
      left_join(
        A_event %>%
          filter(abx_class %in% CLASSES[1:2]) %>%
          select(infant_id, abx_day, abx_class,
                 any_pm7, any_positive_pm7),
        by = c("infant_id", "abx_day", "abx_class")
      ) %>%
      group_by(abx_class, season, z05_flag) %>%
      summarise(
        n_event_days = n(),
        n_tested = sum(any_pm7, na.rm = TRUE),
        n_positive = sum(any_positive_pm7, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        n_negative = n_tested - n_positive,
        n_not_tested = n_event_days - n_tested,
        pct_tested = round(100 * n_tested / n_event_days, 1),
        pct_not_tested = round(100 * n_not_tested / n_event_days, 1)
      ) %>%
      arrange(abx_class, season, desc(z05_flag))

    wcsv(
      z05_events,
      "03_antibiotic_Z051_RSV_testing_combined.csv"
    )

  } else warning('Skipping Z05.1: October diagnosis RDS missing required columns.')
} else invisible(NULL)

wilson <- function(x,n,z=1.96) {
  p=x/n; d=1+z*z/n
  cbind(lo=(p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/d,
        hi=(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/d)
}
plot_dat <- weekly %>% filter(n_tested_infants>0)
if (nrow(plot_dat)) {
  ci <- wilson(plot_dat$n_positive_infants,plot_dat$n_tested_infants)
  plot_dat$lo <- 100*ci[,'lo']; plot_dat$hi <- 100*ci[,'hi']
  p <- ggplot(plot_dat,aes(x=week_of_life,y=positivity_pct,color=season,group=season))+
    geom_line(linewidth=.8)+geom_point(size=2)+
    geom_errorbar(aes(ymin=lo,ymax=hi),width=.12,alpha=.7)+
    scale_x_continuous(breaks=seq_len(MAX_WEEK),labels=paste0('Week ',seq_len(MAX_WEEK)))+
    labs(x='Week of life (week 1 = days 0–6)',y='RSV test positivity (%)',color=NULL,
         title='RSV test positivity by week of life',
         subtitle='Unique infants tested per week; Wilson 95% confidence intervals')+
    theme_bw(base_size=11)+theme(axis.text.x=element_text(angle=35,hjust=1))
  ggsave(file.path(OUT,'04_weekly_RSV_test_positivity.png'),p,width=9,height=5,dpi=300)
}

wcsv(tibble(metric=c('january_cohort_infants','cohort_with_official_new_id',
    'linked_lab_rows_age0_to180','interpretable_test_days_age0_to180',
    'first30_amino_event_days','first30_vanc_event_days',
    'birth_date_column'),value=as.character(c(nrow(cohort),sum(cohort$has_official_new_id),
      nrow(lab),sum(test_days$has_interpretable),sum(events$abx_class==CLASSES[1]),
      sum(events$abx_class==CLASSES[2]),BIRTH_DATE_COLUMN))),
    '00_analysis_QC.csv')
