# Large-scale propensity score (L1-regularized logistic regression on all
# baseline covariates), equipoise check, 1:1 caliper matching, overlap
# weights, and covariate balance diagnostics.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

baseline <- load_derived("baseline")
covariate_labels <- load_derived("covariate_labels")
covariate_ids <- setdiff(names(baseline), c("person_id", "treatment"))

# ---- Propensity model -------------------------------------------------------
x <- scale(as.matrix(baseline[covariate_ids]))
set.seed(2027)
ps_fit <- glmnet::cv.glmnet(x, baseline$treatment, family = "binomial", alpha = 1, nfolds = 10)
baseline$ps <- as.numeric(predict(ps_fit, newx = x, s = "lambda.min", type = "response"))

# Preference score (Walker et al., 2013) rescales the PS so that 0.5 means
# no preference for either drug given the covariates, regardless of market share.
p_target <- mean(baseline$treatment)
baseline$preference <- plogis(qlogis(baseline$ps) - qlogis(p_target))
equipoise <- mean(baseline$preference >= 0.3 & baseline$preference <= 0.7)

selected <- coef(ps_fit, s = "lambda.min")[-1, 1]
ps_coefficients <- tibble(covariate_id = names(selected), coefficient = selected) |>
  filter(coefficient != 0) |>
  left_join(covariate_labels, by = "covariate_id") |>
  arrange(desc(abs(coefficient)))

# ---- 1:1 matching -----------------------------------------------------------
# Each diclofenac user (the smaller cohort) is matched to one celecoxib user
# within 0.2 SD of the logit PS, without replacement. Many patients share a
# PS (most covariates are binary), and MatchIt breaks ties at random.
set.seed(2027)
matching <- MatchIt::matchit(
  treatment ~ age,
  data = baseline,
  distance = qlogis(baseline$ps),
  method = "nearest",
  estimand = "ATC",
  caliper = study$caliper,
  std.caliper = TRUE
)
matched <- MatchIt::match.data(matching, subclass = "pair_id") |>
  select(person_id, treatment, ps, pair_id)

# ---- Overlap weights (sensitivity analysis) ---------------------------------
# Weights each person by the probability of receiving the other drug, which
# targets patients in clinical equipoise and needs no trimming.
baseline$overlap_weight <- WeightIt::get_w_from_ps(baseline$ps, baseline$treatment, estimand = "ATO")

# ---- Balance diagnostics ----------------------------------------------------
is_matched <- baseline$person_id %in% matched$person_id
balance <- covariate_labels |>
  rowwise() |>
  mutate(
    prevalence_or_mean_celecoxib  = mean(baseline[[covariate_id]][baseline$treatment == 1]),
    prevalence_or_mean_diclofenac = mean(baseline[[covariate_id]][baseline$treatment == 0]),
    smd_before   = smd(baseline[[covariate_id]], baseline$treatment),
    smd_matched  = smd(baseline[[covariate_id]][is_matched], baseline$treatment[is_matched]),
    smd_weighted = smd(baseline[[covariate_id]], baseline$treatment, baseline$overlap_weight)
  ) |>
  ungroup() |>
  arrange(desc(abs(smd_before)))

diagnostics <- tibble(
  metric = c(
    "Covariates in PS model", "Covariates selected by LASSO",
    "Share of population in equipoise (preference 0.3-0.7)",
    "Matched pairs", "Celecoxib users unmatched", "Diclofenac users unmatched",
    "Max |SMD| before adjustment", "Max |SMD| after matching", "Max |SMD| after overlap weighting",
    "Covariates with |SMD| > 0.1 after matching"
  ),
  value = c(
    length(covariate_ids), nrow(ps_coefficients),
    round(equipoise, 3),
    nrow(matched) / 2,
    sum(baseline$treatment == 1) - sum(matched$treatment == 1),
    sum(baseline$treatment == 0) - sum(matched$treatment == 0),
    round(max(abs(balance$smd_before)), 3),
    round(max(abs(balance$smd_matched)), 3),
    round(max(abs(balance$smd_weighted)), 3),
    sum(abs(balance$smd_matched) > 0.1)
  )
)

# ---- Table 1 ----------------------------------------------------------------
# Demographics plus the 10 covariates that were most imbalanced before matching.
table1_ids <- union(
  intersect(c("age", "female", "race_white", "race_black", "race_asian", "hispanic"), covariate_ids),
  head(balance$covariate_id, 10)
)
continuous <- c("age", "index_year", "visit_count")
# n (%) for binary covariates, mean (SD) for continuous ones. Suppressed cells
# hide the percentage too, since it would reveal the count.
describe <- function(values, covariate_id) {
  if (covariate_id %in% continuous) return(sprintf("%.1f (%.1f)", mean(values), sd(values)))
  n <- suppress_small_cells(sum(values))
  if (is.na(n)) sprintf("<%d", study$min_cell) else sprintf("%d (%.1f%%)", n, 100 * mean(values))
}
table1 <- balance |>
  filter(covariate_id %in% table1_ids) |>
  rowwise() |>
  mutate(
    celecoxib_before   = describe(baseline[[covariate_id]][baseline$treatment == 1], covariate_id),
    diclofenac_before  = describe(baseline[[covariate_id]][baseline$treatment == 0], covariate_id),
    celecoxib_matched  = describe(baseline[[covariate_id]][is_matched & baseline$treatment == 1], covariate_id),
    diclofenac_matched = describe(baseline[[covariate_id]][is_matched & baseline$treatment == 0], covariate_id)
  ) |>
  ungroup() |>
  transmute(
    characteristic = covariate_name,
    celecoxib_before, diclofenac_before, smd_before = round(smd_before, 3),
    celecoxib_matched, diclofenac_matched, smd_matched = round(smd_matched, 3)
  )

# ---- Figures ----------------------------------------------------------------
arm_label <- function(treatment) factor(ifelse(treatment == 1, "Celecoxib", "Diclofenac"), names(arm_colors))

ps_plot <- ggplot(baseline, aes(preference, fill = arm_label(treatment))) +
  annotate("rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf, alpha = 0.08) +
  geom_density(alpha = 0.55, color = NA) +
  scale_fill_manual(values = arm_colors, name = NULL) +
  labs(x = "Preference score", y = "Density",
       title = "Propensity score overlap",
       subtitle = sprintf("%.0f%% of patients fall in the equipoise region (shaded)", 100 * equipoise)) +
  theme_study()
ggsave(file.path(paths$results, "ps_distribution.png"), ps_plot, width = 6.5, height = 4, dpi = 200)

balance_plot <- ggplot(balance, aes(abs(smd_before), abs(smd_matched))) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, color = "grey80") +
  geom_point(alpha = 0.6, color = arm_colors[["Celecoxib"]]) +
  coord_equal(xlim = c(0, max(abs(balance$smd_before))), ylim = c(0, max(abs(balance$smd_before)))) +
  labs(x = "|SMD| before matching", y = "|SMD| after matching",
       title = "Covariate balance",
       subtitle = sprintf("%d covariates; dashed line marks the 0.1 threshold", nrow(balance))) +
  theme_study()
ggsave(file.path(paths$results, "covariate_balance.png"), balance_plot, width = 5, height = 5, dpi = 200)

# ---- Save -------------------------------------------------------------------
save_derived(matched, "matched")
save_derived(select(baseline, person_id, treatment, ps, overlap_weight), "weighted")
write.csv(diagnostics, file.path(paths$results, "ps_diagnostics.csv"), row.names = FALSE)
write.csv(table1, file.path(paths$results, "table1.csv"), row.names = FALSE)
write.csv(mutate(balance, across(starts_with(c("smd", "prev")), ~ round(.x, 4))),
          file.path(paths$results, "covariate_balance.csv"), row.names = FALSE)
write.csv(mutate(ps_coefficients, coefficient = round(coefficient, 4)),
          file.path(paths$results, "ps_model_coefficients.csv"), row.names = FALSE)
print(diagnostics)
