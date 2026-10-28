# Builds the celecoxib (target) and diclofenac (comparator) new-user cohorts
# and records how many people each eligibility criterion removes.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

con <- connect_omop()
dbExecute(con, render_sql(
  "cohort_candidates.sql",
  target_concept_id     = study$target_concept_id,
  comparator_concept_id = study$comparator_concept_id,
  outcome_concept_id    = study$outcome_concept_id
))
candidates <- dbReadTable(con, "cohort_candidates")

criteria <- list(
  "First exposure to celecoxib or diclofenac"         = quote(TRUE),
  "Index date within an observation period"           = quote(!is.na(prior_obs_days)),
  ">= 365 days of observation before index"           = quote(prior_obs_days >= study$washout_days),
  "Age >= 18 at index"                                = quote(age >= study$min_age),
  "Exposed to only one of the two study drugs"        = quote(!used_both_drugs),
  "No GI bleed on or before index"                    = quote(!prior_outcome)
)

remaining <- candidates
attrition <- vector("list", length(criteria))
for (i in seq_along(criteria)) {
  remaining <- filter(remaining, !!criteria[[i]])
  attrition[[i]] <- tibble(
    step       = i - 1,
    criterion  = names(criteria)[i],
    celecoxib  = sum(remaining$treatment == 1),
    diclofenac = sum(remaining$treatment == 0)
  )
}
attrition <- bind_rows(attrition)

study_population <- remaining |>
  select(person_id, treatment, index_date, obs_end_date, death_date, age)
dbWriteTable(con, "study_population", study_population, overwrite = TRUE)
dbDisconnect(con, shutdown = TRUE)

write.csv(attrition, file.path(paths$results, "attrition.csv"), row.names = FALSE)
print(attrition)
