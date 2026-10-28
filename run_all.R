# Runs the full study end to end. Patient-level data stay in data/ (git-ignored);
# only aggregate tables and figures are written to results/.

scripts <- c(
  "R/01_load_omop.R",
  "R/02_build_cohorts.R",
  "R/03_covariates.R",
  "R/04_propensity_score.R",
  "R/05_outcome_models.R",
  "R/06_negative_controls.R"
)

for (script in scripts) {
  message("\n== ", script, " ==")
  source(here::here(script), local = new.env())
}
