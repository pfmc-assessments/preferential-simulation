library(RTMB)
library(sdmTMB)
library(fmesher)
library(nwfscSurvey)
library(indexwc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(1234)

# Simulation parameters
n_reps <- 25
n_per_year <- 1000 # number of points per year

#############################################################
# Shortspine WCGBTS data
#############################################################
if (!file.exists("simulation_paper/shortspine_wcgbts.rds")) {
  d <- nwfscSurvey::pull_catch(
    survey      = "NWFSC.Combo",
    common_name = "shortspine thornyhead"
  ) |> dplyr::filter(Year %in% 2012:2022)

  names(d) <- tolower(names(d))
  haul <- nwfscSurvey::pull_haul(survey = "NWFSC.Combo")
  haul <- dplyr::filter(haul, trawl_id %in% d$trawl_id)
  d <- dplyr::left_join(d, haul[, c("trawl_id", "area_swept_ha_der")])
  d <- dplyr::select(
    d, year, total_catch_wt_kg, area_swept_ha_der,
    latitude_dd, longitude_dd
  )
  d <- sdmTMB::add_utm_columns(d, ll_names = c("longitude_dd", "latitude_dd"))
  d$fyear <- as.factor(d$year)
  d$log_effort <- log(d$area_swept_ha)
  saveRDS(d, "simulation_paper/shortspine_wcgbts.rds")
}

d <- readRDS("simulation_paper/shortspine_wcgbts.rds")

years <- sort(unique(d$year))
n_years <- length(years)

grid <- sdmTMB::replicate_df(
  indexwc::california_current_grid,
  time_name = "year", time_values = years
)
grid <- dplyr::rename(grid, X = x, Y = y) |>
  dplyr::filter(area_km2_WCGBTS == 4)
grid$fyear <- as.factor(grid$year)
grid_year1 <- dplyr::filter(grid, year == years[1])
grid_year1$cell_id <- seq_len(nrow(grid_year1)) # stable cell index
n_grid <- nrow(grid_year1)

#############################################################
# Baseline sdmTMB (raw data)
#############################################################

# The assessment models use 20km as cutoff, this is intentionally increasing resolution
mesh_baseline <- sdmTMB::make_mesh(d, xy_cols = c("X", "Y"), cutoff = 15)

fit_baseline <- sdmTMB::sdmTMB(
  total_catch_wt_kg ~ 0 + fyear,
  data           = d,
  offset         = log(d$area_swept_ha),
  mesh           = mesh_baseline,
  time           = "year",
  spatiotemporal = list("iid", "iid"),
  family         = delta_gamma()
)

index_baseline <- sdmTMB::get_index(
  predict(fit_baseline, grid, return_tmb_object = TRUE)
)

# calculate the CV from the Gamma, and transform it to log_sd for lognormal
pars <- tidy(fit_baseline, "ran_pars", 2)
gamma_cv <- 1 / sqrt(pars[which(pars$term == "phi"), "estimate"][[1]]) # 1.2
sigma_sim <- sqrt(log(1 + gamma_cv^2))

#############################################################
# Sampling intensity field
#############################################################

set.seed(1234)
samp_locs <- grid_year1[, c("X", "Y")]
samp_locs$year <- years[1]
mesh_samp <- sdmTMB::make_mesh(samp_locs, xy_cols = c("X", "Y"), cutoff = 10)

sim_samp <- sdmTMB::sdmTMB_simulate(
  formula = ~1, data = samp_locs, mesh = mesh_samp,
  family = gaussian(), range = 25, sigma_O = 1.0, sigma_E = 0,
  phi = 0.01, B = 0, seed = 99
)
grid_year1$sampling_score_scaled <- as.numeric(scale(sim_samp$eta))

sim_anom <- sdmTMB::sdmTMB_simulate(
  formula = ~1, data = samp_locs, mesh = mesh_samp,
  family = gaussian(), range = 40, sigma_O = 1.0, sigma_E = 0,
  phi = 0.01, B = 0, seed = 777
)
grid_year1$anom_scaled <- as.numeric(scale(sim_anom$eta))

#############################################################
# Distance field
#############################################################
ports_ll <- data.frame(
  lat = c(
    46.1809, 44.6359, 43.3586, 40.8021, 39.0086, 37.7749,
    36.6002, 35.3734, 34.4308, 33.7490, 32.7157
  ),
  lon = c(
    -123.8008, -124.0562, -124.2158, -124.2720, -123.8061, -122.4194,
    -121.8863, -120.8623, -119.6850, -118.2437, -117.1611
  )
)
ports_utm <- sdmTMB::add_utm_columns(ports_ll,
  ll_names = c("lon", "lat"),
  utm_crs = 32610
) # same
pc <- as.matrix(ports_utm[, c("X", "Y")])
nn <- FNN::get.knnx(pc, as.matrix(grid_year1[, c("X", "Y")]), k = 1)$nn.dist[, 1]
grid_year1$dist_port_scaled <- as.numeric(scale(nn))

# save the static sampling_score field (drawn once)
# This is the observed non-density covariate at each grid cell. It is identical across
# reps and scenarios, so save once.
saveRDS(
  grid_year1[, c("cell_id", "X", "Y", "anom_scaled", "sampling_score_scaled", "dist_port_scaled")],
  "simulation_paper/sampling_grid.rds"
)

#############################################################
# Estimation mesh -- about 580 vertices
#############################################################

mesh <- sdmTMB::make_mesh(grid_year1, xy_cols = c("X", "Y"), cutoff = 20)
spde <- fmesher::fm_fem(mesh$mesh)
n_vertices <- mesh$mesh$n

# Projection from mesh -> every grid cell (built once; grid is fixed)
A_grid <- fmesher::fm_basis(
  mesh$mesh,
  loc = as.matrix(grid_year1[, c("X", "Y")])
)

# Full-grid cell areas for ADREPORT index.
# The full-grid prediction uses grid_year1, so A_grid doubles as A_pred
area_pred_full <- grid_year1$area_km2_WCGBTS

#############################################################
# Ground-truth from sdmTMB model above
#############################################################

pred_grid <- predict(fit_baseline, newdata = grid)
pred_grid$true_density <- plogis(pred_grid$est1) * exp(pred_grid$est2)
pred_grid <- dplyr::left_join(
  pred_grid,
  grid_year1[, c("X", "Y", "cell_id", "sampling_score_scaled", "dist_port_scaled", "anom_scaled")],
  by = c("X", "Y")
)
pred_grid$log_density_scaled <-
  as.numeric(scale(log(pred_grid$true_density + 1e-10)))

# Time-varying preferential sampling:
#   If time_varying_b = TRUE, the preferential coefficient follows a linear
#   trend (0 in year 1 and increases to b_max in the final year) plus AR(1) process
#   deviations. If FALSE, b is CONSTANT at b_max in every year (level-bias
#   regime). b_max indexes the scenario; b_max = 0 is the random control.
time_varying_b <- TRUE # FALSE = constant b

b_max_grid <- c(0, 3, 5) # c(0, 1, 2, 3, 5)
gamma_samp_true <- 1.0 # effect on sampling score field (fixed across all sims)
gamma_dist_true <- 0 #-1.0 # effect of distance on logit of sampling
theta_anom <- 1.0
c0_samp <- -2.0 # arbitrary constant, b_max_grid is interpretable using this intercept
phi_b <- 0.7 # AR(1) parameter
sigma_b <- 0.3 # SD of AR(1) devs

# realized per-year b. Constant regime: b_t = b_max every year. Time-varying
# is a linear ramp mu_t plus stationary AR(1) deviations, lower bound at 0.
# Stochastic in the time-varying case, so each call draws a new trajectory.
b_traj <- function(b_max, years, time_varying = time_varying_b,
                   phi = phi_b, sigma = sigma_b) {
  n <- length(years)
  if (!time_varying) {
    return(setNames(rep(b_max, n), as.character(years))) # constant b
  }
  mu <- if (n == 1) {
    b_max
  } else {
    b_max * (seq_len(n) - 1) / (n - 1) # linear trend
  }
  dev <- numeric(n)
  if (b_max > 0) {
    e <- rnorm(n, 0, sigma) # random devs
    dev[1] <- e[1] / sqrt(1 - phi^2) # stationary init
    for (t in 2:n) dev[t] <- phi * dev[t - 1] + e[t]
  }
  setNames(pmax(mu + dev, 0), as.character(years))
}

# scenario list
scenarios <- setNames(as.list(b_max_grid), paste0("b", b_max_grid))


#############################################################
# Simulate observations (samples grid cells as areal process).
# Sampling weights get calculated using the realized b_t trajectory
# (linear trend + AR(1)); returns the trajectory so each
# replicate b_t can be saved.
#############################################################

simulate_pref_obs <- function(truth, b_max, years, gamma_samp, gamma_dist, c0,
                              n_per_year, sigma_sim, psi, theta, dc = FALSE) {
  b_t <- b_traj(b_max, years) # realized trajectory this rep
  result <- NULL
  for (yr in years) {
    cells <- truth[truth$year == yr, ]
    n_take <- min(n_per_year, nrow(cells))

    if (b_max == 0) {
      w <- rep(1, nrow(cells)) # random-sampling control

      idx <- sample(nrow(cells), size = n_take, replace = TRUE, prob = w)
      s <- cells[idx, ]

      z <- rbinom(n_take, size = 1, prob = plogis(s$est1))
      meanlog <- s$est2 - sigma_sim^2 / 2
      pos <- rlnorm(n_take, meanlog = meanlog, sdlog = sigma_sim)
    } else {
      if (dc == FALSE) {
        # logit_p is affecting which cells are sampled
        # here logit_p is at minimum a function of c0, the sampling intensity
        # field, and log_density_scaled
        logit_p <- c0 + gamma_samp * cells$sampling_score_scaled +
          gamma_dist * cells$dist_port_scaled +
          b_t[as.character(yr)] * cells$log_density_scaled
      } else {
        # DC data-generating model: selection driven by the spatial catch anomaly
        # (anom_scaled) instead of density. The anomaly also inflates the catch
        # (see meanlog below)
        logit_p <- c0 + gamma_samp * cells$sampling_score_scaled +
          gamma_dist * cells$dist_port_scaled +
          b_t[as.character(yr)] * cells$anom_scaled
      }

      # pick cells to sample
      w <- plogis(logit_p)
      idx <- sample(nrow(cells), size = n_take, replace = TRUE, prob = w)
      s <- cells[idx, ]

      # sample occurrence / not in these cells from the density model
      z <- rbinom(n_take, size = 1, prob = plogis(s$est1))
      # positive catch: true density + anomaly (dc only, so the anomaly that
      # drove selection also inflates the catch -> a confounder)
      #anom_effect <- ifelse(dc == TRUE, theta_anom * s$anom_scaled, 0) {
      #  meanlog <- s$est2 + anom_effect - sigma_sim^2 / 2
      #}
      anom_effect <- if (dc) theta_anom * s$anom_scaled else 0
      meanlog <- s$est2 + anom_effect - sigma_sim^2 / 2
      pos <- rlnorm(n_take, meanlog = meanlog, sdlog = sigma_sim)
    }
    result <- rbind(result, data.frame(
      year = yr,
      cell_id = s$cell_id, # stable grid_year1 cell index
      X = s$X,
      Y = s$Y,
      total_catch_wt_kg = z * pos,
      area_swept_ha = 1,
      log_effort = 0,
      sampling_score = s$sampling_score_scaled,
      anom = s$anom_scaled
    ))
  }
  attr(result, "b_t") <- b_t # b_t trajectory
  result
}


for (scenario_name in names(scenarios)) {
  print(scenario_name)
  dat_clean_list <- list()
  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100)

    dat_clean <- simulate_pref_obs(
      pred_grid, scenarios[[scenario_name]], years,
      gamma_samp_true, gamma_dist_true, c0_samp, n_per_year, sigma_sim,
      dc = FALSE
    )
    dat_clean_list[[rep]] <- dat_clean
  }

  saveRDS(dat_clean_list, file = paste0("simulation_paper/data/", scenario_name, ".rds"))
}

# generate alternative datasets from discrete choice approach
for (scenario_name in names(scenarios)) {
  print(scenario_name)
  dat_clean_list <- list()
  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100)

    dat_clean <- simulate_pref_obs(
      pred_grid, scenarios[[scenario_name]], years,
      gamma_samp_true, gamma_dist_true, c0_samp, n_per_year, sigma_sim,
      dc = TRUE
    )
    dat_clean_list[[rep]] <- dat_clean
  }

  saveRDS(dat_clean_list, file = paste0("simulation_paper/data/", scenario_name, "_dc.rds"))
}
