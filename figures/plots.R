library(tidyverse)
library(scales)
library(viridis)
library(dplyr)

scenario_levels <- c("b0", "b3", "b5")
b_end <- c(b0 = 0, b3 = 3, b5 = 5)

# shared styling
my_cols <- viridis::viridis(3, begin = 0.2, end = 0.85, option = "magma")
model_cols <- c(truth = "grey20", naive = my_cols[1], joint = my_cols[2], dc = my_cols[3])
model_labs <- c(truth = "Truth", naive = "Naive", joint = "Joint", dc = "Discrete choice")
model_lty <- c(truth = "dashed", naive = "solid", joint = "solid", dc = "solid")
yr_breaks <- seq(2012, 2022, 4)

# Density driven data
res <- readRDS("simulation_paper/pref_sampling_results_tvb.rds")
resdc <- readRDS("simulation_paper/dchoice_index_tvb.rds") |>
  group_by(scenario, model, rep, grid) |>
  arrange(year, .by_group = TRUE) |>
  mutate(scaled_est = est / est[1]) |>
  ungroup()

results <- res$results
truth_idx <- res$truth

idx1 <- bind_rows(
  results |> filter(model %in% c("naive", "joint")),
  resdc |> filter(model == "dc")
) |>
  filter(grid == "full", converged, scenario %in% scenario_levels) |>
  mutate(scenario = factor(scenario, scenario_levels), dgp = "Density-driven")

index_summary1 <- idx1 |>
  group_by(scenario, model, year) |>
  summarise(
    med = median(scaled_est, na.rm = TRUE),
    lo = quantile(scaled_est, 0.025, na.rm = TRUE),
    hi = quantile(scaled_est, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

truth_facet1 <- crossing(
  scenario = factor(scenario_levels, scenario_levels),
  truth_idx |> transmute(year, med = scaled_est)
) |>
  mutate(model = "truth", lo = NA_real_, hi = NA_real_)

df1 <- bind_rows(index_summary1, truth_facet1) |> mutate(dgp = "Density-driven")

truth_1 <- truth_idx |>
  arrange(year) |>
  transmute(year, truth_scaled = scaled_est, dgp = "Density-driven")

# Anomaly-driven data
res <- readRDS("simulation_paper/pref_sampling_results_tvb_dc.rds")
resdc <- readRDS("simulation_paper/dchoice_index_tvb_dc.rds") |>
  group_by(scenario, model, rep, grid) |>
  arrange(year, .by_group = TRUE) |>
  mutate(scaled_est = est / est[1]) |>
  ungroup()

results <- res$results
truth_idx <- res$truth

idx2 <- bind_rows(
  results |> filter(model %in% c("naive", "joint")),
  resdc |> filter(model == "dc")
) |>
  filter(grid == "full", converged, scenario %in% scenario_levels) |>
  mutate(scenario = factor(scenario, scenario_levels), dgp = "Anomaly-driven")

index_summary2 <- idx2 |>
  group_by(scenario, model, year) |>
  summarise(
    med = median(scaled_est, na.rm = TRUE),
    lo = quantile(scaled_est, 0.025, na.rm = TRUE),
    hi = quantile(scaled_est, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

truth_facet2 <- crossing(
  scenario = factor(scenario_levels, scenario_levels),
  truth_idx |> transmute(year, med = scaled_est)
) |>
  mutate(model = "truth", lo = NA_real_, hi = NA_real_)

df2 <- bind_rows(index_summary2, truth_facet2) |> mutate(dgp = "Anomaly-driven")

truth_2 <- truth_idx |>
  arrange(year) |>
  transmute(year, truth_scaled = scaled_est, dgp = "Anomaly-driven")

# Combine
idx_all <- bind_rows(idx1, idx2)
truth_all <- bind_rows(truth_1, truth_2)

lines_df <- bind_rows(df1, df2) |>
  mutate(
    model = factor(model, levels = c("truth", "naive", "joint", "dc")),
    dgp = factor(dgp, levels = c("Density-driven", "Anomaly-driven")),
    scenario_labs = factor(paste0("b = ", b_end[as.character(scenario)]),
      levels = paste0("b = ", b_end[scenario_levels])
    )
  )

##############################################################################
# FIGURE 3 — Index time series vs truth
##############################################################################
fig1 <- ggplot(lines_df, aes(year, med, colour = model)) +
  geom_ribbon(
    data = filter(lines_df, model != "truth"),
    aes(ymin = lo, ymax = hi, fill = model),
    alpha = 0.15, colour = NA
  ) +
  geom_line(aes(linetype = model), linewidth = 0.8) +
  facet_grid(dgp ~ scenario_labs) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  scale_fill_manual(values = model_cols, guide = "none") +
  scale_linetype_manual(values = model_lty, labels = model_labs, name = NULL) +
  scale_x_continuous(breaks = yr_breaks) +
  labs(x = "Year", y = "Index (scaled to first year)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )

ggsave("simulation_paper/figures/fig3_index_timeseries_combined.png",
  fig1,
  width = 7.5, height = 5
)


##############################################################################
# FIGURE 2 — Index accuracy: log-scale RMSE of the relative index vs truth
##############################################################################
rmse_replicates <- idx_all |>
  left_join(truth_all, by = c("year", "dgp")) |>
  filter(!is.na(scaled_est), !is.na(truth_scaled), scaled_est > 0) |>
  group_by(dgp, scenario, model, rep) |>
  summarise(
    rmse = sqrt(mean((log(scaled_est) - log(truth_scaled))^2)),
    .groups = "drop"
  )

rmse_summary <- rmse_replicates |>
  group_by(dgp, scenario, model) |>
  summarise(
    med = median(rmse),
    lo = quantile(rmse, 0.025),
    hi = quantile(rmse, 0.975),
    .groups = "drop"
  ) |>
  mutate(
    b_end = b_end[as.character(scenario)],
    model = factor(model, c("naive", "joint", "dc")),
    dgp   = factor(dgp, levels = c("Density-driven", "Anomaly-driven"))
  )

rmse_summary <- rmse_summary |>
  transmute(dgp, scenario, model, b_end, med, lo, hi,
    metric = "Index RMSE (log scale)"
  )

truth_slope <- truth_all |>
  group_by(dgp) |>
  summarise(
    true_slope = coef(lm(log(truth_scaled) ~ I(year - min(year))))[2],
    .groups = "drop"
  )

trend_rep <- idx_all |>
  filter(scaled_est > 0) |>
  group_by(dgp, scenario, model, rep) |>
  filter(n() >= 3) |>
  summarise(
    slope = coef(lm(log(scaled_est) ~ I(year - min(year))))[2],
    .groups = "drop"
  ) |>
  left_join(truth_slope, by = "dgp") |>
  mutate(trend_err = slope - true_slope)

trend_summ <- trend_rep |>
  group_by(dgp, scenario, model) |>
  summarise(
    med = median(trend_err),
    lo = quantile(trend_err, 0.025),
    hi = quantile(trend_err, 0.975),
    .groups = "drop"
  ) |>
  mutate(
    b_end = b_end[as.character(scenario)],
    model = factor(model, c("naive", "joint", "dc")),
    dgp = factor(dgp, levels = c("Density-driven", "Anomaly-driven")),
    metric = "Trend error (slope \u2212 truth)"
  )

fig2_df <- bind_rows(rmse_summary, trend_summ) |>
  mutate(metric = factor(metric,
    levels = c(
      "Index RMSE (log scale)",
      "Trend error (slope \u2212 truth)"
    )
  ))

# zero reference line in the trend-error row
zero_line <- data.frame(
  metric = factor("Trend error (slope \u2212 truth)",
    levels = levels(fig2_df$metric)
  ),
  yint = 0
)


fig2 <- ggplot(fig2_df, aes(b_end, med, colour = model)) +
  geom_hline(
    data = zero_line, aes(yintercept = yint),
    linetype = "dotted", colour = "grey45"
  ) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
    position = position_dodge(width = 0.3),
    linewidth = 0.7, fatten = 3
  ) +
  facet_grid(metric ~ dgp, scales = "free_y", switch = "y") +
  scale_colour_manual(
    values = model_cols[c("naive", "joint", "dc")],
    labels = model_labs[c("naive", "joint", "dc")], name = NULL
  ) +
  scale_x_continuous(breaks = b_end) +
  labs(x = "Preferential strength scenario (b)", y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    legend.position = "bottom"
  )

ggsave("simulation_paper/figures/fig4_index_rmse.png",
  fig2,
  width = 7, height = 6
)


##############################################################################
# FIGURE 4 — Indices in stock assessment
##############################################################################

get_estimate <- function(ss3_out, out_table, label) {
  purrr::map_dbl(ss3_out, \(x) {
    filter(x[[out_table]], Label == label) |>
      pull(Value)
  })
}
quad_results_nested_ss3 <- readRDS(
  "simulation_paper/quadrature_ss3_results.rds"
) |>
  dplyr::filter(scenario %in% c("b0", "b3", "b5"))

fig4 <- quad_results_nested_ss3 |>
  mutate(
    Recr_unfished = get_estimate(ss3_out, "derived_quants", "Recr_unfished"),
    #  logr0 = get_estimate(ss3_out, 'parameters', 'SR_LN(R0)'),
    ofl_2025 = get_estimate(ss3_out, "derived_quants", "OFLCatch_2025"),
    ofl_2034 = get_estimate(ss3_out, "derived_quants", "OFLCatch_2034"),
    recdev_2008 = get_estimate(ss3_out, "parameters", "Main_RecrDev_2008"),
    scenario = as.numeric(stringr::str_remove(scenario, "b")),
    model = factor(model, c("naive", "joint"))
  ) |>
  select(-ss3_tab, -ss3_out) |>
  tidyr::pivot_longer(
    Recr_unfished:recdev_2008,
    names_to = "metric",
    values_to = "estimate"
  ) |>
  mutate(
    metric = factor(
      case_when(
        grepl("2025", metric) ~ "OFL[2025]",
        grepl("2034", metric) ~ "OFL[2034]",
        grepl("Recr", metric) ~ "R[0]",
        grepl("recdev", metric) ~ "epsilon[2008]"
      ),
      levels = c("R[0]", "epsilon[2008]", "OFL[2025]", "OFL[2034]")
    )
  ) |>
  group_by(scenario, model, metric) |>
  summarise(
    med = median(estimate),
    lo.mid = quantile(estimate, 0.25),
    hi.mid = quantile(estimate, 0.75),
    lo = quantile(estimate, 0.025),
    hi = quantile(estimate, 0.975),
    .groups = "drop_last"
  ) |>
  # fig4 <- fig4_tab |>
  ggplot(aes(x = scenario, col = model)) +
  geom_hline(
    data = \(df) {
      df |>
        filter(scenario == 0) |>
        select(metric, ref = med)
    },
    aes(yintercept = ref),
    col = "red"
  ) +
  geom_pointrange(
    aes(y = med, ymin = lo, ymax = hi),
    position = position_dodge(width = 0.75),
    linewidth = 0.7,
    fatten = 3
  ) +
  geom_linerange(
    aes(ymin = lo.mid, ymax = hi.mid),
    position = position_dodge(width = 0.75),
    linewidth = 1.5,
  ) +
  facet_wrap(~metric, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(
    values = model_cols[c("naive", "joint")],
    labels = model_labs[c("naive", "joint")]
  ) +
  scale_x_continuous(breaks = b_end, name = "Scenario") +
  theme_bw() +
  theme(
    # legend.position = 'bottom',
    strip.background = element_rect(fill = "white")
  ) +
  ylab("Estimate")

ggsave(
  filename = "simulation_paper/figures/fig4_ss3_results.png",
  fig4,
  width = 5,
  height = 4
)
