# Negative control outcomes: re-runs the matched analysis for conditions that
# neither NSAID should cause. Their estimates should center on HR = 1; any
# systematic shift measures residual bias, which empirical calibration then
# folds into the GI bleeding estimate's p-value and confidence interval.
#
# Negative controls are too rare within 90 days, so both they and the GI
# bleeding estimate being calibrated use the intention-to-treat window. In
# this dataset every GI bleed occurs within 90 days of index, so the ITT and
# primary estimates coincide.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

matched <- load_derived("matched")

con <- connect_omop(read_only = TRUE)
primary <- get_follow_up(con, study$outcome_concept_id, study$tar_itt) |>
  inner_join(matched, by = "person_id") |>
  fit_cox("matched")
concept_names <- dbGetQuery(con, sprintf(
  "SELECT concept_id, concept_name FROM concept WHERE concept_id IN (%s)",
  paste(study$negative_controls, collapse = ", ")
))

nc_estimates <- lapply(study$negative_controls, function(concept_id) {
  data <- get_follow_up(con, concept_id, study$tar_itt) |>
    filter(!prior_outcome) |>
    inner_join(matched, by = "person_id")
  # Too few events give unstable estimates that would distort the null fit.
  if (sum(data$event) < 5 || length(unique(data$treatment[data$event == 1])) < 2) return(NULL)
  mutate(fit_cox(data, "matched"), concept_id = concept_id, .before = 1)
})
dbDisconnect(con, shutdown = TRUE)

nc_estimates <- bind_rows(nc_estimates) |>
  left_join(concept_names, by = "concept_id") |>
  relocate(concept_name, .after = concept_id)

# ---- Empirical calibration --------------------------------------------------
null <- EmpiricalCalibration::fitNull(nc_estimates$log_hr, nc_estimates$se_log_hr)
error_model <- EmpiricalCalibration::fitSystematicErrorModel(
  nc_estimates$log_hr, nc_estimates$se_log_hr, trueLogRr = rep(0, nrow(nc_estimates))
)
calibrated_p <- EmpiricalCalibration::calibrateP(null, primary$log_hr, primary$se_log_hr)
calibrated_ci <- EmpiricalCalibration::calibrateConfidenceInterval(
  primary$log_hr, primary$se_log_hr, error_model
)

calibration <- tibble(
  estimate = c("Uncalibrated", "Calibrated"),
  hr       = c(primary$hr, exp(calibrated_ci$logRr)),
  ci_lower = c(primary$ci_lower, exp(calibrated_ci$logLb95Rr)),
  ci_upper = c(primary$ci_upper, exp(calibrated_ci$logUb95Rr)),
  p        = c(primary$p, calibrated_p)
)
null_distribution <- tibble(
  parameter = c("mean (log HR)", "sd (log HR)", "negative controls used"),
  value     = c(null[["mean"]], null[["sd"]], nrow(nc_estimates))
)

EmpiricalCalibration::plotCalibrationEffect(
  logRrNegatives   = nc_estimates$log_hr,
  seLogRrNegatives = nc_estimates$se_log_hr,
  logRrPositives   = primary$log_hr,
  seLogRrPositives = primary$se_log_hr,
  null             = null,
  xLabel           = "Hazard ratio",
  title            = "Negative controls (blue) and the GI bleeding estimate (yellow)",
  fileName         = file.path(paths$results, "negative_controls.png")
)

write.csv(round_estimates(nc_estimates),
          file.path(paths$results, "negative_control_estimates.csv"), row.names = FALSE)
write.csv(mutate(calibration, across(where(is.double), ~ signif(.x, 4))),
          file.path(paths$results, "calibrated_estimate.csv"), row.names = FALSE)
write.csv(mutate(null_distribution, value = signif(value, 4)),
          file.path(paths$results, "systematic_error.csv"), row.names = FALSE)
print(select(nc_estimates, concept_name, hr, ci_lower, ci_upper))
print(null_distribution)
print(calibration)
