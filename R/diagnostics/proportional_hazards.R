

library(tidyverse)
library(survival)

needed <- c("model_data_primary", "model_data_fullyear", "fit_ag")
missing_objects <- needed[!vapply(needed, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects/functions: ",
    paste(missing_objects, collapse = ", "),
    "\nRun AG_Cox_postindex_antibiotic_calendar_month_strata(4).R first."
  )
}

stopifnot("all_cause" %in% names(model_data_primary))
stopifnot("all_cause" %in% names(model_data_fullyear))

if (exists("cox_primary") && "all_cause" %in% names(cox_primary)) {
  fit_rsv_180 <- cox_primary$all_cause
} else {
  fit_rsv_180 <- fit_ag(model_data_primary$all_cause)
}

fit_full_180 <- fit_ag(model_data_fullyear$all_cause)

restrict_age_window <- function(dat, lower, upper) {

  dat %>%
    filter(tstop > lower, tstart < upper) %>%
    mutate(

      event = as.integer(event == 1 & tstop <= upper & tstop > lower),
      tstart = pmax(tstart, lower),
      tstop  = pmin(tstop, upper)
    ) %>%
    filter(tstop > tstart)
}

rsv_1_90    <- restrict_age_window(model_data_primary$all_cause, 0, 90)
rsv_91_180  <- restrict_age_window(model_data_primary$all_cause, 90, 180)

full_1_90   <- restrict_age_window(model_data_fullyear$all_cause, 0, 90)
full_91_180 <- restrict_age_window(model_data_fullyear$all_cause, 90, 180)

fit_rsv_1_90    <- fit_ag(rsv_1_90)
fit_rsv_91_180  <- fit_ag(rsv_91_180)

fit_full_1_90   <- fit_ag(full_1_90)
fit_full_91_180 <- fit_ag(full_91_180)

fits <- list(
  "RSV seasons: 1-180 days" = fit_rsv_180,
  "RSV seasons: 1-90 days" = fit_rsv_1_90,
  "RSV seasons: 91-180 days" = fit_rsv_91_180,
  "Full year: 1-180 days" = fit_full_180,
  "Full year: 1-90 days" = fit_full_1_90,
  "Full year: 91-180 days" = fit_full_91_180
)

zph_tests <- purrr::imap(fits, function(fit, analysis) {
  z <- cox.zph(fit, transform = "km", terms = FALSE)

  as.data.frame(z$table) %>%
    tibble::rownames_to_column("term") %>%
    mutate(analysis = analysis, .before = 1)
})

schoenfeld_tests <- bind_rows(zph_tests) %>%
  mutate(
    p_formatted = case_when(
      p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p)
    )
  )
print(schoenfeld_tests)

write.csv(
  schoenfeld_tests,
  "schoenfeld_tests_allcause.csv",
  row.names = FALSE
)

exposure_tests <- schoenfeld_tests %>%
  filter(
    str_detect(str_to_lower(term), "maternal|nirsevimab") |
      term == "GLOBAL"
  )
print(exposure_tests)

write.csv(
  exposure_tests,
  "schoenfeld_tests_exposure_only.csv",
  row.names = FALSE
)

get_exposure_indices <- function(z) {
  nm <- colnames(z$y)
  maternal <- grep("maternal", nm, ignore.case = TRUE)
  nirsevimab <- grep("nirsevimab", nm, ignore.case = TRUE)
  if (length(maternal) != 1L || length(nirsevimab) != 1L) {
    print(nm)
    stop("Could not uniquely identify maternal and nirsevimab coefficients.")
  }
  c(maternal = maternal, nirsevimab = nirsevimab)
}

plot_fits <- list(
  fit_rsv_180, fit_rsv_180, fit_full_180, fit_full_180,
  fit_rsv_1_90, fit_rsv_1_90, fit_full_1_90, fit_full_1_90,
  fit_rsv_91_180, fit_rsv_91_180, fit_full_91_180, fit_full_91_180
)

zph_plot <- lapply(plot_fits, function(fit) {
  cox.zph(fit, transform = "identity", terms = FALSE)
})

panel_exposure <- rep(c("maternal", "nirsevimab"), 6)
panel_letters <- LETTERS[1:12]
panel_context <- rep(
  c("RSV seasons", "RSV seasons", "Full year", "Full year"), 3
)
panel_windows <- rep(c("days 1-180", "days 1-90", "days 91-180"), each = 4)
panel_exposure_label <- ifelse(
  panel_exposure == "maternal", "Maternal RSVpreF", "Nirsevimab"
)

get_row_ylim <- function(panel_ids) {
  values <- unlist(lapply(panel_ids, function(i) {
    z <- zph_plot[[i]]
    j <- unname(get_exposure_indices(z)[panel_exposure[i]])
    z$y[, j]
  }), use.names = FALSE)
  values <- values[is.finite(values)]
  if (!length(values)) stop("No finite Schoenfeld residuals in row.")
  lim <- range(values)
  pad <- max(diff(lim) * 0.06, 0.5)
  lim + c(-pad, pad)
}
row_ylims <- lapply(0:2, function(r) get_row_ylim((4*r + 1):(4*r + 4)))

output_pdf <- "Figure_SX_Schoenfeld_diagnostics_combined_v8.pdf"
pdf(output_pdf, width = 7, height = 8.8 * (2/3), pointsize = 8,
    family = "Helvetica", useDingbats = FALSE)

layout(matrix(1:12, nrow = 3, byrow = TRUE),
       widths = rep(1, 4), heights = rep(1, 3))

par(lwd = 0.65, tck = -0.025, mgp = c(1.7, 0.35, 0),
    cex.axis = 0.74)

for (i in seq_len(12)) {
  row_id <- (i - 1) %/% 4 + 1
  col_id <- (i - 1) %% 4 + 1
  z <- zph_plot[[i]]
  j <- unname(get_exposure_indices(z)[panel_exposure[i]])
  ylim <- row_ylims[[row_id]]

  par(mar = c(if (row_id == 3) 2.3 else 1.35,
              if (col_id == 1) 3.15 else 0.65, 3.2, 0.35))

  plot(z, var = j, resid = FALSE, se = TRUE,
       axes = FALSE, ann = FALSE, xlab = "", ylab = "",
       xlim = range(z$x[is.finite(z$x)]), ylim = ylim,
       col = "black")

  valid <- is.finite(z$x) & is.finite(z$y[, j])
  points(z$x[valid], z$y[valid, j],
         pch = 1, cex = 0.35,
         col = grDevices::adjustcolor("black", alpha.f = 0.30))

  x_at <- pretty(range(z$x[is.finite(z$x)]), n = 4)
  x_at <- x_at[x_at >= min(z$x, na.rm = TRUE) &
               x_at <= max(z$x, na.rm = TRUE)]
  y_at <- pretty(ylim, n = 5)
  y_at <- y_at[y_at >= ylim[1] & y_at <= ylim[2]]

  axis(1, at = x_at, labels = x_at, lwd = 0,
       lwd.ticks = 0.65, cex.axis = 0.72)
  axis(2, at = y_at,
       labels = if (col_id == 1) y_at else FALSE,
       las = 1, lwd = 0, lwd.ticks = 0.65,
       mgp = c(1.7, 0.85, 0), cex.axis = 0.72)
  box(bty = "l", lwd = 0.65)

  if (row_id == 3) {
    mtext("Infant age (days)", side = 1, line = 1.3, cex = 0.73)
  }
  if (col_id == 1) {
    mtext("Coefficient", side = 2, line = 1.7, cex = 0.74)
  }

  header <- sprintf("%s. %s, %s\n%s", panel_letters[i],
                    panel_context[i], panel_windows[i],
                    panel_exposure_label[i])
  mtext(header, side = 3, line = 1.0, cex = 0.74,
        font = 2, adj = 0)
}

dev.off()
model_checks <- tibble(
  analysis = c(
    "RSV seasons: 1-180 days",
    "RSV seasons: 1-90 days",
    "RSV seasons: 91-180 days",
    "Full year: 1-180 days",
    "Full year: 1-90 days",
    "Full year: 91-180 days"
  ),
  n_rows = c(
    nrow(model_data_primary$all_cause),
    nrow(rsv_1_90),
    nrow(rsv_91_180),
    nrow(model_data_fullyear$all_cause),
    nrow(full_1_90),
    nrow(full_91_180)
  ),
  infants = c(
    n_distinct(model_data_primary$all_cause$infant_study_id),
    n_distinct(rsv_1_90$infant_study_id),
    n_distinct(rsv_91_180$infant_study_id),
    n_distinct(model_data_fullyear$all_cause$infant_study_id),
    n_distinct(full_1_90$infant_study_id),
    n_distinct(full_91_180$infant_study_id)
  ),
  events = c(
    sum(model_data_primary$all_cause$event),
    sum(rsv_1_90$event),
    sum(rsv_91_180$event),
    sum(model_data_fullyear$all_cause$event),
    sum(full_1_90$event),
    sum(full_91_180$event)
  )
)
print(model_checks)

write.csv(
  model_checks,
  "schoenfeld_model_checks.csv",
  row.names = FALSE
)
