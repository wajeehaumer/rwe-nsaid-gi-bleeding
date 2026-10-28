suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
  library(ggplot2)
})

paths <- list(
  db      = here::here("data", "omop.duckdb"),
  derived = here::here("data", "derived"),
  results = here::here("results")
)
dir.create(paths$derived, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)

connect_omop <- function(read_only = FALSE) {
  dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = paths$db, read_only = read_only)
}

# Minimal SqlRender-style templating: replaces @name with the supplied value.
render_sql <- function(file, ...) {
  sql <- paste(readLines(here::here("sql", file), warn = FALSE), collapse = "\n")
  params <- list(...)
  for (name in names(params)) {
    sql <- gsub(paste0("@", name, "\\b"), paste(params[[name]], collapse = ", "), sql)
  }
  sql
}

# Mask nonzero counts below the threshold most data use agreements require
# before aggregate results leave the secure environment.
suppress_small_cells <- function(x, min_cell = study$min_cell) {
  ifelse(x > 0 & x < min_cell, NA, x)
}

# Standardized mean difference for a binary or continuous covariate.
smd <- function(x, treatment, weights = rep(1, length(x))) {
  w1 <- weights[treatment == 1]; x1 <- x[treatment == 1]
  w0 <- weights[treatment == 0]; x0 <- x[treatment == 0]
  m1 <- weighted.mean(x1, w1); m0 <- weighted.mean(x0, w0)
  v1 <- sum(w1 * (x1 - m1)^2) / sum(w1)
  v0 <- sum(w0 * (x0 - m0)^2) / sum(w0)
  pooled_sd <- sqrt((v1 + v0) / 2)
  if (pooled_sd == 0) return(0)
  (m1 - m0) / pooled_sd
}

# Time-to-first-outcome data for one outcome concept and time-at-risk window.
# People whose follow-up ends before the window opens are dropped; people
# with the outcome on or before index are flagged so callers can exclude them.
get_follow_up <- function(con, outcome_concept_id, tar) {
  dbGetQuery(con, render_sql(
    "time_at_risk.sql",
    outcome_concept_id = outcome_concept_id,
    tar_start = tar[1],
    tar_end   = tar[2]
  )) |>
    filter(tar_end >= tar_start) |>
    mutate(
      event = as.integer(!is.na(outcome_date)),
      time  = as.numeric(coalesce(outcome_date, tar_end) - index_date)
    )
}

tar_label <- function(tar) {
  if (tar[2] >= 36500) sprintf("Day %d to end of observation", tar[1]) else sprintf("Days %d-%d", tar[1], tar[2])
}

# Cox proportional hazards estimate of celecoxib vs diclofenac with robust
# standard errors: clustered on matched pair, or sandwich for weights.
fit_cox <- function(data, analysis = c("crude", "matched", "overlap")) {
  analysis <- match.arg(analysis)
  fit <- switch(analysis,
    crude   = survival::coxph(survival::Surv(time, event) ~ treatment, data = data),
    matched = survival::coxph(survival::Surv(time, event) ~ treatment, data = data, cluster = pair_id),
    overlap = survival::coxph(survival::Surv(time, event) ~ treatment, data = data,
                              weights = overlap_weight, robust = TRUE)
  )
  log_hr <- unname(coef(fit)["treatment"])
  se <- sqrt(fit$var[1, 1])
  arm <- function(t) filter(data, treatment == t)
  person_years <- function(d) sum(d$time) / 365.25
  tibble(
    analysis,
    persons_celecoxib  = suppress_small_cells(nrow(arm(1))),
    persons_diclofenac = suppress_small_cells(nrow(arm(0))),
    events_celecoxib   = suppress_small_cells(sum(arm(1)$event)),
    events_diclofenac  = suppress_small_cells(sum(arm(0)$event)),
    ir_per_1000py_celecoxib  = 1000 * sum(arm(1)$event) / person_years(arm(1)),
    ir_per_1000py_diclofenac = 1000 * sum(arm(0)$event) / person_years(arm(0)),
    hr = exp(log_hr),
    ci_lower = exp(log_hr - qnorm(0.975) * se),
    ci_upper = exp(log_hr + qnorm(0.975) * se),
    log_hr,
    se_log_hr = se,
    p = 2 * pnorm(-abs(log_hr / se))
  )
}

round_estimates <- function(estimates) {
  mutate(estimates, across(c(starts_with("ir_"), hr, ci_lower, ci_upper, log_hr, se_log_hr, p), ~ signif(.x, 4)))
}

save_derived <- function(x, name) saveRDS(x, file.path(paths$derived, paste0(name, ".rds")))
load_derived <- function(name) readRDS(file.path(paths$derived, paste0(name, ".rds")))

theme_study <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

arm_colors <- c(Celecoxib = "#2a6f97", Diclofenac = "#c8553d")
