

library(tidyverse)

required_objects <- c(
  "primary_antibiotic_records_eligible",
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

required_cols <- c(
  "infant_study_id",
  "abx_days",
  "abx_class",
  diagnosis_flags
)

missing_cols <- setdiff(
  required_cols,
  names(primary_antibiotic_records_eligible)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in primary_antibiotic_records_eligible: ",
    paste(missing_cols, collapse = ", "),
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
  all(diagnosis_flags %in% names(DX_LABELS))
)

AMINO_CLASS <- "Aminoglycoside Antibiotics"
VANC_CLASS  <- "Vancomycin Antibiotics And Derivatives"

target_classes <- c(
  AMINO_CLASS,
  VANC_CLASS
)

observed_classes <- unique(
  primary_antibiotic_records_eligible$abx_class
)

missing_classes <- setdiff(
  target_classes,
  observed_classes
)

if (length(missing_classes) > 0) {
  stop(
    "Target antibiotic classes not found: ",
    paste(missing_classes, collapse = ", "),
    call. = FALSE
  )
}

c6_records_1_30 <- primary_antibiotic_records_eligible %>%
  filter(
    !is.na(abx_days),
    abx_days >= 1,
    abx_days <= 30
  )

class_event_keys <- c6_records_1_30 %>%
  filter(
    !is.na(abx_class),
    abx_class %in% target_classes
  ) %>%
  transmute(
    infant_study_id,
    abx_days,
    antibiotic_class = case_when(
      abx_class == AMINO_CLASS ~ "Aminoglycoside",
      abx_class == VANC_CLASS  ~ "Vancomycin"
    )
  ) %>%
  distinct()

diagnosis_by_event_day <- c6_records_1_30 %>%
  group_by(
    infant_study_id,
    abx_days
  ) %>%
  summarise(
    across(
      all_of(diagnosis_flags),
      ~ as.integer(any(.x == 1, na.rm = TRUE))
    ),
    .groups = "drop"
  )

c6_events <- class_event_keys %>%
  left_join(
    diagnosis_by_event_day,
    by = c(
      "infant_study_id",
      "abx_days"
    )
  )

stopifnot(
  nrow(c6_events) ==
    nrow(
      distinct(
        c6_events,
        infant_study_id,
        abx_days,
        antibiotic_class
      )
    )
)

c6_event_counts <- c6_events %>%
  group_by(
    antibiotic_class
  ) %>%
  summarise(
    n_events = n(),
    n_infants = n_distinct(infant_study_id),
    .groups = "drop"
  )
print(c6_event_counts)

summarize_c6_diagnoses <- function(data, class_name) {

  d <- data %>%
    filter(
      antibiotic_class == class_name
    )

  denominator <- nrow(d)

  if (denominator == 0) {
    stop(
      "No events found for ",
      class_name,
      " during days 1-30.",
      call. = FALSE
    )
  }

  dx_rows <- d %>%
    summarise(
      across(
        all_of(diagnosis_flags),
        ~ sum(.x == 1, na.rm = TRUE)
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "diagnosis_flag",
      values_to = "n"
    ) %>%
    mutate(
      diagnosis = unname(
        DX_LABELS[diagnosis_flag]
      ),
      denominator = denominator,
      percent = 100 * n / denominator
    ) %>%
    select(
      diagnosis_flag,
      diagnosis,
      n,
      denominator,
      percent
    )

  no_dx <- d %>%
    rowwise() %>%
    mutate(
      n_positive_dx = sum(
        c_across(
          all_of(diagnosis_flags)
        ) == 1,
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    summarise(
      n = sum(n_positive_dx == 0)
    ) %>%
    mutate(
      diagnosis_flag = "none_prespecified",
      diagnosis =
        "No prespecified coded diagnosis category identified",
      denominator = denominator,
      percent = 100 * n / denominator
    ) %>%
    select(
      diagnosis_flag,
      diagnosis,
      n,
      denominator,
      percent
    )

  bind_rows(
    dx_rows,
    no_dx
  ) %>%
    mutate(
      antibiotic_class = class_name,
      `n (%)` = sprintf(
        "%s (%.1f%%)",
        scales::comma(n),
        percent
      )
    ) %>%
    select(
      antibiotic_class,
      diagnosis_flag,
      diagnosis,
      n,
      denominator,
      percent,
      `n (%)`
    )
}

c6_long <- bind_rows(
  summarize_c6_diagnoses(
    c6_events,
    "Aminoglycoside"
  ),
  summarize_c6_diagnoses(
    c6_events,
    "Vancomycin"
  )
)

diagnosis_order <- c6_long %>%
  filter(
    diagnosis_flag != "none_prespecified"
  ) %>%
  group_by(
    diagnosis_flag,
    diagnosis
  ) %>%
  summarise(
    total_n = sum(n),
    .groups = "drop"
  ) %>%
  arrange(
    desc(total_n)
  ) %>%
  pull(
    diagnosis_flag
  )

diagnosis_order <- c(
  diagnosis_order,
  "none_prespecified"
)

c6_table <- c6_long %>%
  select(
    diagnosis_flag,
    diagnosis,
    antibiotic_class,
    `n (%)`
  ) %>%
  pivot_wider(
    names_from = antibiotic_class,
    values_from = `n (%)`
  ) %>%
  mutate(
    diagnosis_flag = factor(
      diagnosis_flag,
      levels = diagnosis_order
    )
  ) %>%
  arrange(
    diagnosis_flag
  ) %>%
  select(
    -diagnosis_flag
  )
print(
  c6_table,
  n = Inf,
  width = Inf
)

c6_dx_burden <- c6_events %>%
  rowwise() %>%
  mutate(
    n_prespecified_dx = sum(
      c_across(
        all_of(diagnosis_flags)
      ) == 1,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  count(
    antibiotic_class,
    n_prespecified_dx,
    name = "n_events"
  ) %>%
  group_by(
    antibiotic_class
  ) %>%
  mutate(
    percent = 100 * n_events / sum(n_events)
  ) %>%
  ungroup()
print(
  c6_dx_burden,
  n = Inf
)

write.csv(
  c6_table,
  "C6_amino_vanc_diagnoses_days1_30.csv",
  row.names = FALSE
)

write.csv(
  c6_long,
  "C6_amino_vanc_diagnoses_days1_30_long.csv",
  row.names = FALSE
)

write.csv(
  c6_event_counts,
  "C6_amino_vanc_event_counts_days1_30.csv",
  row.names = FALSE
)

write.csv(
  c6_dx_burden,
  "C6_amino_vanc_diagnosis_burden_days1_30.csv",
  row.names = FALSE
)

invisible(NULL)
