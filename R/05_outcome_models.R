# Hazard ratios for GI bleeding, celecoxib vs diclofenac: crude, PS-matched
# (primary), and overlap-weighted, across the primary and sensitivity
# time-at-risk windows. Also the Kaplan-Meier curve and E-value.

source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))

matched  <- load_derived("matched")
weighted <- load_derived("weighted")

con <- connect_omop(read_only = TRUE)
windows <- c(list(study$tar_primary), study$tar_sensitivity, list(study$tar_itt))

estimates <- lapply(windows, function(tar) {
  follow_up <- get_follow_up(con, study$outcome_concept_id, tar)
  label <- tar_label(tar)
  bind_rows(
    fit_cox(inner_join(follow_up, weighted, by = "person_id"), "crude"),
    fit_cox(inner_join(follow_up, matched,  by = "person_id"), "matched"),
    fit_cox(inner_join(follow_up, weighted, by = "person_id"), "overlap")
  ) |>
    mutate(time_at_risk = label, primary_window = identical(tar, study$tar_primary), .before = 1)
})
estimates <- bind_rows(estimates)

primary_follow_up <- inner_join(get_follow_up(con, study$outcome_concept_id, study$tar_primary),
                                matched, by = "person_id")
dbDisconnect(con, shutdown = TRUE)

# E-value: how strong an unmeasured confounder would need to be, on the risk
# ratio scale with both treatment and outcome, to explain away the estimate.
# The outcome is not rare in this population, so rare = FALSE.
primary <- filter(estimates, primary_window, analysis == "matched")
e_values <- EValue::evalues.HR(primary$hr, primary$ci_lower, primary$ci_upper, rare = FALSE)

# ---- Kaplan-Meier curve (matched cohort, primary window) --------------------
km <- survival::survfit(survival::Surv(time, event) ~ treatment, data = primary_follow_up)
km_data <- broom::tidy(km) |>
  mutate(arm = factor(ifelse(strata == "treatment=1", "Celecoxib", "Diclofenac"), names(arm_colors)))

km_plot <- ggplot(km_data, aes(time, 1 - estimate, color = arm, fill = arm)) +
  geom_ribbon(aes(ymin = 1 - conf.high, ymax = 1 - conf.low), alpha = 0.15, color = NA, stat = "identity") +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = arm_colors, name = NULL) +
  scale_fill_manual(values = arm_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days since treatment initiation", y = "Cumulative incidence of GI bleeding",
       title = "GI bleeding after NSAID initiation (PS-matched cohort)",
       subtitle = sprintf("Matched HR %.2f (95%% CI %.2f-%.2f)", primary$hr, primary$ci_lower, primary$ci_upper)) +
  theme_study()
ggsave(file.path(paths$results, "km_cumulative_incidence.png"), km_plot, width = 6.5, height = 4.5, dpi = 200)

# ---- Forest plot across analyses --------------------------------------------
forest_data <- estimates |>
  mutate(
    analysis = factor(analysis, c("crude", "matched", "overlap"),
                      c("Crude", "PS-matched 1:1", "Overlap-weighted")),
    time_at_risk = factor(time_at_risk, unique(time_at_risk))
  )
forest_plot <- ggplot(forest_data, aes(hr, analysis, xmin = ci_lower, xmax = ci_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(color = analysis == "PS-matched 1:1"), show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = arm_colors[["Celecoxib"]], `FALSE` = "grey35")) +
  scale_x_log10() +
  facet_wrap(~ time_at_risk, ncol = 1) +
  labs(x = "Hazard ratio, celecoxib vs diclofenac (log scale)", y = NULL,
       title = "GI bleeding: estimates across analyses and time-at-risk windows") +
  theme_study()
ggsave(file.path(paths$results, "forest_plot.png"), forest_plot, width = 6.5, height = 6, dpi = 200)

write.csv(round_estimates(estimates),
          file.path(paths$results, "hazard_ratios.csv"), row.names = FALSE)
write.csv(as.data.frame(e_values), file.path(paths$results, "e_value.csv"))
print(select(estimates, time_at_risk, analysis, hr, ci_lower, ci_upper))
print(e_values)
