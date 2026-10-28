# Assembles the baseline covariate matrix: demographics plus data-driven
# condition, drug, and utilization covariates (in the spirit of OHDSI
# FeatureExtraction), keeping those above a minimum prevalence.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

con <- connect_omop()
dbExecute(con, render_sql(
  "covariates.sql",
  window_start          = study$covariate_window[1],
  window_end            = study$covariate_window[2],
  target_concept_id     = study$target_concept_id,
  comparator_concept_id = study$comparator_concept_id
))
covariates_long <- dbReadTable(con, "covariates_long")

demographics <- dbGetQuery(con, "
  SELECT sp.person_id,
         sp.age,
         CAST(p.gender_concept_id = 8532 AS DOUBLE)        AS female,
         CAST(p.race_concept_id = 8527 AS DOUBLE)          AS race_white,
         CAST(p.race_concept_id = 8516 AS DOUBLE)          AS race_black,
         CAST(p.race_concept_id = 8515 AS DOUBLE)          AS race_asian,
         CAST(p.ethnicity_concept_id = 38003563 AS DOUBLE) AS hispanic,
         year(sp.index_date)                               AS index_year
  FROM study_population sp
  JOIN person p ON p.person_id = sp.person_id")
population <- dbReadTable(con, "study_population")
dbDisconnect(con, shutdown = TRUE)

# Keep binary covariates seen in at least min_covariate_prevalence of the
# population; rarer ones add noise to the propensity model and can't be
# balanced meaningfully.
n <- nrow(population)
keep <- covariates_long |>
  count(covariate_id, covariate_name) |>
  filter(covariate_id == "visit_count" | n / !!n >= study$min_covariate_prevalence)

wide <- covariates_long |>
  semi_join(keep, by = "covariate_id") |>
  select(person_id, covariate_id, value) |>
  tidyr::pivot_wider(names_from = covariate_id, values_from = value, values_fill = 0)
# Fixed column order too, so the LASSO fit doesn't depend on query row order.
wide <- wide[c("person_id", sort(setdiff(names(wide), "person_id")))]

baseline <- population |>
  select(person_id, treatment) |>
  left_join(demographics, by = "person_id") |>
  left_join(wide, by = "person_id") |>
  mutate(across(-c(person_id, treatment), ~ coalesce(.x, 0))) |>
  # Database row order isn't guaranteed; a fixed order keeps cross-validation
  # folds and matching reproducible from run to run.
  arrange(person_id)

# Drop covariates with no variation (e.g., a diagnosis everyone has).
constant <- names(baseline)[vapply(baseline, function(x) length(unique(x)) == 1, logical(1))]
baseline <- select(baseline, -all_of(constant))

covariate_labels <- bind_rows(
  tibble(
    covariate_id   = c("age", "female", "race_white", "race_black", "race_asian", "hispanic", "index_year"),
    covariate_name = c("Age (years)", "Female", "Race: White", "Race: Black", "Race: Asian", "Hispanic", "Index year")
  ),
  distinct(keep, covariate_id, covariate_name)
) |>
  filter(covariate_id %in% names(baseline))

save_derived(baseline, "baseline")
save_derived(covariate_labels, "covariate_labels")
message(nrow(baseline), " persons x ", ncol(baseline) - 2, " covariates",
        if (length(constant)) paste0(" (dropped constant: ", paste(constant, collapse = ", "), ")"))
