

library(tidyverse)

required_objects <- c(
  "primary_event_days_with_diagnoses_eligible",
  "diagnosis_flags"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Run the main AG analysis and updated supplement script first. Missing: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

DX_LABELS <- c(
  asthma_allergy_yn      = "Asthma/allergy",
  bronchitis_yn          = "Bronchitis/bronchiolitis",
  fever_yn               = "Fever",
  influenza_yn           = "Influenza",
  bacterial_infection_yn = "Bacterial infection (misc.)",
  other_infection_yn     = "Other infection (misc.)",
  sup_om_yn              = "Acute/suppurative otitis media",
  nsup_om                = "Non-suppurative otitis media",
  pharyngitis_yn         = "Pharyngitis",
  pneumonia_yn           = "Pneumonia",
  sinusitis_yn           = "Sinusitis",
  UTI_yn                 = "Urinary tract infection",
  skin_cond_yn           = "Skin/cutaneous/mucosal condition",
  URI_yn                 = "Viral upper respiratory infection"
)

stopifnot(
  all(diagnosis_flags %in% names(DX_LABELS)),
  all(diagnosis_flags %in%
        names(primary_event_days_with_diagnoses_eligible))
)

c5_event_days <- primary_event_days_with_diagnoses_eligible

stopifnot(
  nrow(c5_event_days) ==
    nrow(distinct(c5_event_days, infant_study_id, abx_days))
)

n_all <- nrow(c5_event_days)
n_uri <- sum(c5_event_days$URI_yn == 1, na.rm = TRUE)
n_bronchitis <- sum(c5_event_days$bronchitis_yn == 1, na.rm = TRUE)
expected_counts <- c(
  all = 4358L,
  uri = 452L,
  bronchitis = 244L
)

observed_counts <- c(
  all = n_all,
  uri = n_uri,
  bronchitis = n_bronchitis
)

if (!all(observed_counts == expected_counts)) {
  warning(
    paste0(
      "C5 validation counts do not match the current supplement. ",
      "Expected all/URI/bronchitis = ",
      paste(expected_counts, collapse = "/"),
      "; observed = ",
      paste(observed_counts, collapse = "/"),
      ". Check that the current main analysis and supplement scripts were ",
      "run in the same R session before interpreting C5 results."
    )
  )
}

diagnosis_flag_missingness <- c5_event_days %>%
  summarise(
    across(
      all_of(diagnosis_flags),
      ~ sum(is.na(.x))
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "diagnosis_flag",
    values_to = "n_missing"
  ) %>%
  mutate(
    diagnosis = unname(DX_LABELS[diagnosis_flag]),
    pct_missing = 100 * n_missing / n_all
  ) %>%
  select(diagnosis_flag, diagnosis, n_missing, pct_missing)
print(diagnosis_flag_missingness, n = Inf)

overall_diagnosis_linkage <- c5_event_days %>%
  rowwise() %>%
  mutate(
    n_positive_dx = sum(c_across(all_of(diagnosis_flags)) == 1, na.rm = TRUE),
    n_missing_dx  = sum(is.na(c_across(all_of(diagnosis_flags)))),
    diagnosis_status = case_when(
      n_positive_dx > 0 ~ "At least one coded diagnosis category identified",
      n_positive_dx == 0 & n_missing_dx == 0 ~
        "No coded diagnosis category identified",
      TRUE ~ "No positive category, but >=1 diagnosis flag missing"
    )
  ) %>%
  ungroup() %>%
  count(diagnosis_status, name = "n") %>%
  mutate(
    denominator = n_all,
    percent = 100 * n / denominator,
    `n (%)` = sprintf("%s (%.1f%%)", scales::comma(n), percent)
  )
print(overall_diagnosis_linkage, n = Inf)

make_cooccurrence_table <- function(data, index_flag, index_label) {

  if (!index_flag %in% diagnosis_flags) {
    stop("Index flag is not in diagnosis_flags: ", index_flag, call. = FALSE)
  }

  index_days <- data %>%
    filter(.data[[index_flag]] == 1)

  denominator <- nrow(index_days)
  other_flags <- setdiff(diagnosis_flags, index_flag)

  if (denominator == 0) {
    stop("No antibiotic-use days found for ", index_label, ".", call. = FALSE)
  }

  co_dx <- map_dfr(
    other_flags,
    function(flag) {
      tibble(
        index_syndrome = index_label,
        index_flag = index_flag,
        denominator = denominator,
        cooccurring_diagnosis = unname(DX_LABELS[flag]),
        cooccurring_flag = flag,
        n = sum(index_days[[flag]] == 1, na.rm = TRUE),
        n_missing_flag = sum(is.na(index_days[[flag]]))
      )
    }
  ) %>%
    mutate(
      percent = 100 * n / denominator,
      `n (%)` = sprintf("%s (%.1f%%)", scales::comma(n), percent)
    )

  no_other <- index_days %>%
    rowwise() %>%
    mutate(
      n_other_positive = sum(
        c_across(all_of(other_flags)) == 1,
        na.rm = TRUE
      ),
      n_other_missing = sum(
        is.na(c_across(all_of(other_flags)))
      )
    ) %>%
    ungroup()

  no_other_complete <- no_other %>%
    filter(n_other_positive == 0, n_other_missing == 0) %>%
    nrow()

  unknown_other <- no_other %>%
    filter(n_other_positive == 0, n_other_missing > 0) %>%
    nrow()

  summary_rows <- tibble(
    index_syndrome = index_label,
    index_flag = index_flag,
    denominator = denominator,
    cooccurring_diagnosis = c(
      "No other coded diagnosis category identified",
      "No other positive category, but >=1 other diagnosis flag missing"
    ),
    cooccurring_flag = c(NA_character_, NA_character_),
    n = c(no_other_complete, unknown_other),
    n_missing_flag = c(NA_integer_, NA_integer_)
  ) %>%
    mutate(
      percent = 100 * n / denominator,
      `n (%)` = sprintf("%s (%.1f%%)", scales::comma(n), percent)
    )

  if (unknown_other == 0) {
    summary_rows <- summary_rows %>%
      filter(
        cooccurring_diagnosis !=
          "No other positive category, but >=1 other diagnosis flag missing"
      )
  }

  bind_rows(co_dx, summary_rows)
}

uri_cooccurrence <- make_cooccurrence_table(
  data = c5_event_days,
  index_flag = "URI_yn",
  index_label = "Viral upper respiratory infection"
)

bronchitis_cooccurrence <- make_cooccurrence_table(
  data = c5_event_days,
  index_flag = "bronchitis_yn",
  index_label = "Bronchitis/bronchiolitis"
)

c5_cooccurrence_long <- bind_rows(
  uri_cooccurrence,
  bronchitis_cooccurrence
)
print(
  uri_cooccurrence %>%
    select(
      cooccurring_diagnosis,
      n,
      denominator,
      percent,
      `n (%)`
    ),
  n = Inf
)
print(
  bronchitis_cooccurrence %>%
    select(
      cooccurring_diagnosis,
      n,
      denominator,
      percent,
      `n (%)`
    ),
  n = Inf
)

display_order <- c(
  "Acute/suppurative otitis media",
  "Non-suppurative otitis media",
  "Pneumonia",
  "Pharyngitis",
  "Sinusitis",
  "Influenza",
  "Urinary tract infection",
  "Skin/cutaneous/mucosal condition",
  "Bacterial infection (misc.)",
  "Other infection (misc.)",
  "Fever",
  "Asthma/allergy",
  "Bronchitis/bronchiolitis",
  "Viral upper respiratory infection",
  "No other coded diagnosis category identified",
  "No other positive category, but >=1 other diagnosis flag missing"
)

c5_table <- c5_cooccurrence_long %>%
  select(
    index_syndrome,
    cooccurring_diagnosis,
    `n (%)`
  ) %>%
  pivot_wider(
    names_from = index_syndrome,
    values_from = `n (%)`
  ) %>%

  complete(
    cooccurring_diagnosis = display_order
  ) %>%
  mutate(
    `Viral upper respiratory infection` = if_else(
      cooccurring_diagnosis == "Viral upper respiratory infection",
      "",
      `Viral upper respiratory infection`
    ),
    `Bronchitis/bronchiolitis` = if_else(
      cooccurring_diagnosis == "Bronchitis/bronchiolitis",
      "",
      `Bronchitis/bronchiolitis`
    ),
    display_order_num = match(cooccurring_diagnosis, display_order)
  ) %>%
  filter(

    !(
      cooccurring_diagnosis ==
        "No other positive category, but >=1 other diagnosis flag missing" &
      is.na(`Viral upper respiratory infection`) &
      is.na(`Bronchitis/bronchiolitis`)
    )
  ) %>%
  arrange(display_order_num) %>%
  select(
    `Co-occurring diagnosis` = cooccurring_diagnosis,
    `Viral upper respiratory infection`,
    `Bronchitis/bronchiolitis`
  )
print(c5_table, n = Inf, width = Inf)
write_csv(
  c5_cooccurrence_long,
  "C5_cooccurring_diagnoses_long.csv"
)

write_csv(
  c5_table,
  "C5_cooccurring_diagnoses_table.csv"
)

write_csv(
  overall_diagnosis_linkage,
  "C5_overall_diagnosis_linkage.csv"
)

write_csv(
  diagnosis_flag_missingness,
  "C5_diagnosis_flag_missingness.csv"
)
