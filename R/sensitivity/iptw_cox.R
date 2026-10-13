

stopifnot(exists('model_data_primary'), exists('model_data_fullyear'),
          exists('cohort_primary'), exists('COVARIATE_TERMS'))
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(data.table)
  library(splines); library(broom)
})
OUT <- 'AG_IPTW_allcause_output_followup_corrected'
dir.create(OUT, showWarnings=FALSE)
EPS <- 0.005
TRUNC <- c(.01,.99)

full <- as.data.table(model_data_fullyear$all_cause)
ins <- as.data.table(model_data_primary$all_cause)
stopifnot(all(c('infant_study_id','mom_study_id','tstart','tstop',
                'event','exposure_arm','calendar_month','index_date') %in% names(full)))

base <- full[order(infant_study_id,tstart), .SD[1], by=infant_study_id]
base[, maternal := as.integer(as.character(exposure_arm)=='a_maternal_vax_only')]

base[, maternal := as.integer(as.character(exposure_maternal)=='maternal_vax')]
base[, infant_id := as.character(infant_study_id)]

dose <- as.data.table(cohort_primary)[,.(infant_study_id, mab_history, mab_days)]
base <- merge(base,dose,by='infant_study_id',all.x=TRUE,sort=FALSE)
base[, dose_day := fifelse(!is.na(mab_history)&mab_history &
                            is.finite(mab_days)&mab_days>=0 & mab_days<180,
                          as.integer(floor(mab_days)), NA_integer_)]

stopifnot(!any(base$maternal==1L & !is.na(base$dose_day)))

baseline_terms <- paste('GENDER + race_cat_infant + serious_newborn_comorbidity +',
 'apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +',
 'median_income_cat + OUTPT + INPT + ED + flu_prior + tdap_prior +',
 'tot_covid + gbs_prior + lb_prior + ns(as.numeric(index_date), df=5)')
ps_m <- glm(as.formula(paste('maternal ~',baseline_terms)),
            data=base, family=binomial())
base[, pm := pmin(1-EPS,pmax(EPS,as.numeric(predict(ps_m,type='response'))))]
p_m <- mean(base$maternal)
base[, w_m := fifelse(maternal==1L,p_m/pm,(1-p_m)/(1-pm))]

nonmat <- base[maternal==0L, .(infant_study_id, infant_id, dose_day,index_date,
 GENDER,race_cat_infant,serious_newborn_comorbidity,apgar_1min,age_index,
 charlson_Wt,bmi_prior,smoking,median_income_cat,OUTPT,INPT,ED,
 flu_prior,tdap_prior,tot_covid,gbs_prior,lb_prior)]

observed <- full[, .(last_stop=max(tstop)), by=infant_study_id]
nonmat <- merge(nonmat,observed,by='infant_study_id',all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(nonmat$last_stop))
nonmat[, last_eligible_day := as.integer(ceiling(last_stop)-1L)]
nonmat[, day_end := pmin(179L,last_eligible_day,
                        fifelse(is.na(dose_day),179L,dose_day))]
days <- nonmat[day_end>=0L][, .(day=0L:day_end), by=names(nonmat)]

observed_intervals <- unique(full[,.(infant_study_id,tstart,tstop)])
setkey(observed_intervals,infant_study_id,tstart,tstop)
day_windows <- unique(days[,.(infant_study_id,day)])
day_windows[, `:=`(day_start=as.numeric(day),day_stop=as.numeric(day+1L))]
setkey(day_windows,infant_study_id,day_start,day_stop)
valid_days <- foverlaps(day_windows,observed_intervals,
                       by.x=c('infant_study_id','day_start','day_stop'),
                       by.y=c('infant_study_id','tstart','tstop'),
                       type='any',nomatch=0L)[pmin(day_stop,tstop)>pmax(day_start,tstart),
                                               .(infant_study_id,day)]
valid_days <- unique(valid_days)
setkey(valid_days,infant_study_id,day)
setkey(days,infant_study_id,day)
days <- days[valid_days,nomatch=0L]
stopifnot(!any(days$day >= days$last_stop))

days[, initiated := as.integer(!is.na(dose_day) & day==dose_day)]
days[, calendar_day := as.numeric(index_date + day)]

utd <- unique(full[as.character(exposure_arm)=='d_neither',
  .(infant_study_id,tstart,utd_status)])
setorder(utd,infant_study_id,tstart)
setkey(utd,infant_study_id,tstart)
days[, lookup_day := as.numeric(day)]
days[, utd_status:=utd[days,on=.(infant_study_id,tstart=lookup_day),roll=TRUE,x.utd_status]]
days[is.na(utd_status),utd_status:='Up to date']
days[, utd_status:=factor(utd_status)]

common <- 'ns(day,df=5) + ns(calendar_day,df=8)'
confounders <- paste('GENDER + race_cat_infant + serious_newborn_comorbidity +',
 'apgar_1min + age_index + charlson_Wt + bmi_prior + smoking +',
 'median_income_cat + OUTPT + INPT + ED + flu_prior + tdap_prior +',
 'tot_covid + gbs_prior + lb_prior + utd_status')
fit_den <- glm(as.formula(paste('initiated ~',common,'+',confounders)),
               data=days,family=binomial())
fit_num <- glm(as.formula(paste('initiated ~',common)),
               data=days,family=binomial())
days[, pd := pmin(1-EPS,pmax(EPS,as.numeric(predict(fit_den,type='response'))))]
days[, pn := pmin(1-EPS,pmax(EPS,as.numeric(predict(fit_num,type='response'))))]

days[, log_increment := fifelse(initiated==1L,log(pn)-log(pd),
                                log1p(-pn)-log1p(-pd))]
setorder(days,infant_study_id,day)
days[, log_cum := cumsum(log_increment),by=infant_study_id]
days[, w_n := exp(pmin(20,pmax(-20,log_cum)))]

lookup_intervals <- unique(full[as.character(exposure_maternal)!='maternal_vax',
                            .(infant_study_id,tstart,tstop)])
lookup_intervals[, start_day:=as.integer(floor(tstart))]
lookup_intervals[, end_day:=as.integer(ceiling(tstop)-1L)]
all_days <- unique(lookup_intervals[end_day>=start_day,
                   .(day=seq.int(start_day,end_day)),
                   by=.(infant_study_id,tstart,tstop)][,.(infant_study_id,day)])
setkey(days,infant_study_id,day)
all_days <- days[all_days,on=.(infant_study_id,day),roll=TRUE,
                 .(infant_study_id=i.infant_study_id,day=i.day,w_n)]
all_days[is.na(w_n),w_n:=1]

expand_ag <- function(md) {
  x <- as.data.table(copy(md)); x[, original_row:=.I]
  x[, start_day:=as.integer(floor(tstart))]
  x[, end_day:=as.integer(ceiling(tstop)-1L)]
  x <- x[end_day>=start_day]
  x <- x[, .(day=seq.int(start_day,end_day)),by=names(x)]
  x[, new_start:=pmax(tstart,day)]
  x[, new_stop:=pmin(tstop,day+1)]
  x <- x[new_stop>new_start]
  x[, event:=as.integer(event==1L & abs(new_stop-tstop)<1e-8)]
  x[, `:=`(tstart=new_start,tstop=new_stop)]
  x[,c('new_start','new_stop'):=NULL]
  x
}

run_period <- function(md,label) {
  x <- expand_ag(md)
  x <- merge(x,base[,.(infant_study_id,maternal,w_m)],by='infant_study_id',all.x=TRUE,sort=FALSE)
  x <- merge(x,all_days,by=c('infant_study_id','day'),all.x=TRUE,sort=FALSE)
  x[maternal==1L,w_n:=1]
  if (anyNA(x$w_m) || anyNA(x$w_n)) stop('Missing IPTW after daily merge')
  x[, sw_raw:=w_m*w_n]

  unique_w <- unique(x[,.(infant_study_id,day,sw_raw)])
  q <- quantile(unique_w$sw_raw,TRUNC,na.rm=TRUE)
  x[, sw:=pmax(q[1],pmin(q[2],sw_raw))]
  x[, exposure_arm:=relevel(droplevels(factor(exposure_arm)),ref='d_neither')]
  f <- Surv(tstart,tstop,event) ~ exposure_arm + strata(calendar_month)
  fit <- coxph(f,data=x,weights=sw,cluster=mom_study_id,
               ties='efron',robust=TRUE)
  res <- broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE)
  res <- res[grepl('^exposure_arm',res$term),]
  res <- transform(res,period=label,HR=estimate,CI_low=conf.low,
                   CI_high=conf.high,VE_pct=100*(1-estimate),
                   VE_CI_low=100*(1-conf.high),VE_CI_high=100*(1-conf.low))
  diagnostics <- data.frame(period=label,rows=nrow(x),
    infants=uniqueN(x$infant_study_id),events=sum(x$event),
    weight_p01=unname(q[1]),weight_p99=unname(q[2]),
    raw_min=min(unique_w$sw_raw),raw_max=max(unique_w$sw_raw),
    truncated_max=max(x$sw),n_unique_infant_days=nrow(unique_w))
  list(results=res,diagnostics=diagnostics)
}

r_season <- run_period(ins,'RSV_season')
r_full <- run_period(full,'Full_year')
r_out <- if (exists('model_data_out_of_season')) run_period(as.data.table(model_data_out_of_season$all_cause),'Non_RSV_season') else NULL
results <- rbindlist(lapply(Filter(Negate(is.null),list(r_season,r_full,r_out)), `[[`, 'results'),fill=TRUE)
diags <- rbindlist(lapply(Filter(Negate(is.null),list(r_season,r_full,r_out)), `[[`, 'diagnostics'),fill=TRUE)
write.csv(results,file.path(OUT,'IPTW_AG_allcause_estimates.csv'),row.names=FALSE)
write.csv(diags,file.path(OUT,'IPTW_weight_diagnostics.csv'),row.names=FALSE)
invisible(NULL); print(results[,c('period','term','HR','CI_low','CI_high','VE_pct','VE_CI_low','VE_CI_high')])
invisible(NULL);print(diags)
