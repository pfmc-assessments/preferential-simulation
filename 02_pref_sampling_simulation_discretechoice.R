library(RTMB)
library(sdmTMB)
library(fmesher)
library(indexwc)
library(FNN)
library(dplyr)

time_varying_b <- TRUE
b_tag <- if (time_varying_b) "tvb" else "constb"
K_nearby <- 25 # candidates per choice set
# pool_size <- 100 # draw the 50 candidates randomly from this many cells (large cells break this, and
# make the code very very slow)
set.seed(1)

# rebuild the same fixed grid / mesh / projections as the sim
scenario_names <- c("b3", "b5") # c("b1", "b2", "b3", "b5")

# read data for first iteration (b1)
init_data <- readRDS(paste0("simulation_paper/data/", "b3", ".rds"))[[1]]
years <- sort(unique(init_data$year))
n_years <- length(years)

# make the prediction grid from the WCGBTS trawl survey grid
grid <- sdmTMB::replicate_df(indexwc::california_current_grid,
  time_name = "year", time_values = years
)
# filter out partial cells, or cells part of other surveys
grid <- dplyr::rename(grid, X = x, Y = y) |>
  dplyr::filter(area_km2_WCGBTS == 4)
grid_year1 <- dplyr::filter(grid, year == years[1])
grid_year1$cell_id <- seq_len(nrow(grid_year1))
n_grid <- nrow(grid_year1)

# make the mesh using sdmTMB, and cutoff distance of 20km
mesh <- sdmTMB::make_mesh(grid_year1, xy_cols = c("X", "Y"), cutoff = 20)
spde <- fmesher::fm_fem(mesh$mesh)
n_vertices <- mesh$mesh$n
A_pred <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(grid_year1[, c("X", "Y")]))
area_pred_full <- grid_year1$area_km2_WCGBTS

# For every grid cell, draw the candidate set from all cells belonging to its
# nearest `n_nodes` vertices (not its nearest cells). This widens the
# choice-set footprint to a few mesh nodes
n_nodes <- 4 # number of nearest mesh vertices whose cells form the pool
grid_coords <- as.matrix(grid_year1[, c("X", "Y")])
mesh_coords <- mesh$mesh$loc[, 1:2]

# which mesh vertex is each grid cell nearest to
cell_node <- FNN::get.knnx(mesh_coords, grid_coords, k = 1)$nn.index[, 1] # n_grid
# cells grouped by mesh vertex
cells_by_node <- split(seq_len(n_grid), cell_node)
# for each cell, its nearest n_nodes mesh vertices
cell_nn_nodes <- FNN::get.knnx(mesh_coords, grid_coords, k = n_nodes)$nn.index # n_grid x n_nodes

nearby_all <- matrix(0L, n_grid, K_nearby)
nearby_all[, 1] <- seq_len(n_grid) # self (chosen cell), always column 1
for (c in seq_len(n_grid)) {
  # pool = all cells belonging to this cell's nearest n_nodes vertices
  pool <- unlist(cells_by_node[as.character(cell_nn_nodes[c, ])], use.names = FALSE)
  pool <- setdiff(pool, c) # exclude self (it's column 1)
  if (length(pool) >= K_nearby - 1) {
    nearby_all[c, 2:K_nearby] <- sample(pool, K_nearby - 1)
  } else {
    # small pool: sample with replacement to fill (rare, edge vertices)
    nearby_all[c, 2:K_nearby] <- sample(pool, K_nearby - 1, replace = TRUE)
  }
}


#############################################################
# RTMB data list
#############################################################
make_data <- function(dat_clean, mesh, spde, A_pred, area_pred, nearby_all) {
  years_v <- sort(unique(dat_clean$year))
  year_idx <- match(dat_clean$year, years_v)
  A_obs <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(dat_clean[, c("X", "Y")]))
  n_grid <- nrow(A_pred)

  # per-observation 50-cell choice set, expanded from the table for each cell by
  # cell_id (chosen cell = column 1). Stored 0-indexed, as in the real-data model.
  nearby_cells <- nearby_all[dat_clean$cell_id, , drop = FALSE] - 1L # n_obs x 50

  list(
    observations = dat_clean$total_catch_wt_kg,
    log_effort = dat_clean$log_effort,
    year_indices = year_idx,
    pos_idx = which(dat_clean$total_catch_wt_kg > 0),
    nearby_cells = nearby_cells, # n_obs x 50 cell indices, 0-indexed (chosen = col 1)
    K_nearby = ncol(nearby_cells),
    n_obs = nrow(dat_clean),
    n_years = length(years_v),
    spde = spde,
    A_obs = A_obs,
    A_pred = A_pred,
    n_grid = n_grid,
    n_vertices = mesh$mesh$n,
    area_pred = area_pred,
    # Currently I'm fixing the range for the lambda SVC fields. I set this so
    # the correction can vary at the same spatial
    # scale as density but cannot drift to a bad range.
    kappa_lambda_fixed = sqrt(8) / 25, # range is 25 km
    # latent eta field: patchy, fixed range (~40 km) matching the true anomaly
    kappa_eta_fixed = sqrt(8) / 40
  )
}

#############################################################
# RTMB parameter list
#############################################################
make_pars <- function(dat_clean, n_vertices, n_years) {
  list(
    year_effects_pos = rep(
      log(mean(dat_clean$total_catch_wt_kg[dat_clean$total_catch_wt_kg > 0])) -
        mean(dat_clean$log_effort),
      n_years
    ),
    year_effects_bin = rep(0, n_years),
    log_sigma_O_pos = log(0.3),
    log_sigma_E_pos = log(0.3),
    log_sigma_O_bin = log(0.3),
    log_sigma_E_bin = log(0.3),
    log_sigma_obs = log(0.5),
    log_kappa = log(sqrt(8) / 300),
    omega_s_pos = rep(0, n_vertices),
    omega_s_bin = rep(0, n_vertices),
    epsilon_st_pos = matrix(0, n_vertices, n_years),
    epsilon_st_bin = matrix(0, n_vertices, n_years),
    log_sigma_alpha = log(0.3),
    alpha_catch = rep(0, n_years),
    b_eta = 0,
    log_sigma_eta = log(0.5),
    omega_eta = rep(0, n_vertices),
    log_sigma_lambda1 = log(0.5),
    log_sigma_lambda2 = log(0.5),
    log_sigma_lambdac = log(0.5),
    omega_lambda1 = rep(0, n_vertices),
    omega_lambda2 = rep(0, n_vertices),
    omega_lambdac = rep(0, n_vertices)
  )
}

#############################################################
# NLL - main function
#############################################################
make_nll_spatiotemporal <- function(tmb_data) {
  function(par) {
    getAll(par, tmb_data)
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    # these could be a single number but splitting it out as I was tuning
    eps_prob <- 1e-12
    eps_dens <- 1e-12
    eps_scale <- 1e-12

    kappa <- exp(log_kappa)
    range_shared <- sqrt(8) / kappa
    sigma_O_pos <- exp(log_sigma_O_pos) + eps_scale
    sigma_E_pos <- exp(log_sigma_E_pos) + eps_scale
    sigma_O_bin <- exp(log_sigma_O_bin) + eps_scale
    sigma_E_bin <- exp(log_sigma_E_bin) + eps_scale
    sigma_obs <- exp(log_sigma_obs) + eps_scale

    Q_base <- kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2
    tau_O_pos <- 1 / (sqrt(4 * pi) * kappa * sigma_O_pos)
    tau_E_pos <- 1 / (sqrt(4 * pi) * kappa * sigma_E_pos)
    tau_O_bin <- 1 / (sqrt(4 * pi) * kappa * sigma_O_bin)
    tau_E_bin <- 1 / (sqrt(4 * pi) * kappa * sigma_E_bin)
    Q_O_pos <- Q_base * tau_O_pos^2
    Q_E_pos <- Q_base * tau_E_pos^2
    Q_O_bin <- Q_base * tau_O_bin^2
    Q_E_bin <- Q_base * tau_E_bin^2

    omega_s_pos %~% dgmrf(0, Q_O_pos)
    omega_s_bin %~% dgmrf(0, Q_O_bin)
    for (t in 1:n_years) epsilon_st_pos[, t] %~% dgmrf(0, Q_E_pos)
    for (t in 1:n_years) epsilon_st_bin[, t] %~% dgmrf(0, Q_E_bin)

    kappa_lambda <- kappa_lambda_fixed
    sigma_lambda1 <- exp(log_sigma_lambda1) + eps_scale
    sigma_lambda2 <- exp(log_sigma_lambda2) + eps_scale
    sigma_lambdac <- exp(log_sigma_lambdac) + eps_scale
    Q_base_lambda <- kappa_lambda^4 * spde$c0 + 2 * kappa_lambda^2 * spde$g1 + spde$g2
    tau_lambda1 <- 1 / (sqrt(4 * pi) * kappa_lambda * sigma_lambda1)
    tau_lambda2 <- 1 / (sqrt(4 * pi) * kappa_lambda * sigma_lambda2)
    tau_lambdac <- 1 / (sqrt(4 * pi) * kappa_lambda * sigma_lambdac)
    omega_lambda1 %~% dgmrf(0, Q_base_lambda * tau_lambda1^2)
    omega_lambda2 %~% dgmrf(0, Q_base_lambda * tau_lambda2^2)
    omega_lambdac %~% dgmrf(0, Q_base_lambda * tau_lambdac^2)

    # latent eta field: patchy (fixed range), shared across years. Enters
    # the utility (so p_chosen is affected) and the positive-catch predictor (so
    # catch residual is affected). Not included in the index, so eta
    # absorbing the anomaly-driven catch inflation is also a correction?
    kappa_eta <- kappa_eta_fixed
    sigma_eta <- exp(log_sigma_eta) + eps_scale
    Q_base_eta <- kappa_eta^4 * spde$c0 + 2 * kappa_eta^2 * spde$g1 + spde$g2
    tau_eta <- 1 / (sqrt(4 * pi) * kappa_eta * sigma_eta)
    omega_eta %~% dgmrf(0, Q_base_eta * tau_eta^2)
    eta_grid <- as.vector(A_pred %*% omega_eta)
    eta_obs <- as.vector(A_obs %*% omega_eta)

    # expected density surface: computed 1x per year (passed to the utility)
    omega_grid_pos <- as.vector(A_pred %*% omega_s_pos)
    omega_grid_bin <- as.vector(A_pred %*% omega_s_bin)
    eps_grid_pos <- A_pred %*% epsilon_st_pos
    eps_grid_bin <- A_pred %*% epsilon_st_bin

    log_expected_density_pred <- matrix(0, n_grid, n_years)
    for (t in 1:n_years) {
      le <- year_effects_bin[t] + omega_grid_bin + eps_grid_bin[, t]
      ld <- year_effects_pos[t] + omega_grid_pos + eps_grid_pos[, t]
      # this is keeping the encounter prob away from 0 so log() can't hit -Inf
      pe <- eps_prob + (1 - 2 * eps_prob) * plogis(le)
      log_expected_density_pred[, t] <- log(pe) + ld
    }

    # year-specific preferential strength via random walk (like b_t in the Conn joint model)
    sigma_alpha <- exp(log_sigma_alpha) + eps_scale
    alpha_catch[1] %~% dnorm(0, 2.0, log = TRUE) # diffuse year-1 prior
    for (t in 2:n_years) {
      alpha_catch[t] %~% dnorm(alpha_catch[t - 1], sigma_alpha, log = TRUE)
    }

    # discrete-choice likelihood over each observation cells 50-cell set
    # utility = alpha * log(expected density) + eta (latent selection field)
    alpha_obs <- alpha_catch[year_indices]
    eVk <- vector("list", K_nearby)
    for (k in 1:K_nearby) {
      nearby_idx <- nearby_cells[, k] + 1
      Vk <- alpha_obs *
        log_expected_density_pred[cbind(nearby_idx, year_indices)] +
        eta_grid[nearby_idx]
      eVk[[k]] <- exp(Vk)
    }
    denom <- eVk[[1]]
    for (k in 2:K_nearby) denom <- denom + eVk[[k]]
    p_chosen <- eVk[[1]] / (denom + eps_dens) # chosen cell = column 1

    # keep p_chosen strictly in (0,1) for dbinom
    p_chosen <- eps_prob + (1 - 2 * eps_prob) * p_chosen

    # choice likelihood
    as.numeric(rep(1, n_obs)) %~% dbinom(size = 1, prob = p_chosen)

    # lambda coefficients vary in space (projected from the SVC fields to each
    # observation): (lambda_c(s) + lambda_1(s)*p + lambda_2(s)*p^2).
    # Allen: this is basically what the economics papers are doing (this imposes a smooth constraint)
    lambdac_i <- as.vector(A_obs %*% omega_lambdac)
    lambda1_i <- as.vector(A_obs %*% omega_lambda1)
    lambda2_i <- as.vector(A_obs %*% omega_lambda2)
    poly_correction <- lambdac_i + lambda1_i * p_chosen + lambda2_i * p_chosen^2

    # process model (delta-lognormal) at the hauls and vectorized. Encounter
    # (binomial) is uncorrected; positive (lognormal) uses the polynomial correction, like Chen et al.
    omega_proj_pos <- as.vector(A_obs %*% omega_s_pos)
    omega_proj_bin <- as.vector(A_obs %*% omega_s_bin)
    eps_proj_pos <- (A_obs %*% epsilon_st_pos)[cbind(seq_len(n_obs), year_indices)]
    eps_proj_bin <- (A_obs %*% epsilon_st_bin)[cbind(seq_len(n_obs), year_indices)]

    logit_encounter <- year_effects_bin[year_indices] +
      omega_proj_bin + eps_proj_bin
    log_abundance_pos <- year_effects_pos[year_indices] +
      omega_proj_pos + eps_proj_pos + eta_obs + log_effort
    log_abundance_corrected <- log_abundance_pos + poly_correction

    REPORT(p_chosen)
    REPORT(log_abundance_pos)
    p_encounter <- plogis(logit_encounter)
    p_encounter <- eps_prob + (1 - 2 * eps_prob) * p_encounter
    as.numeric(observations > 0) %~% dbinom(size = 1, prob = p_encounter)
    log(observations[pos_idx]) %~%
      dnorm(log_abundance_corrected[pos_idx], sigma_obs)

    b_eta %~% dnorm(0, 2.0, log = TRUE)
    log_sigma_eta %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_lambda1 %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_lambda2 %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_lambdac %~% dnorm(log(0.5), 0.5, log = TRUE)
    ADREPORT(sigma_alpha)
    ADREPORT(alpha_catch) # lenght is n_years path
    ADREPORT(b_eta)
    ADREPORT(sigma_eta)
    ADREPORT(sigma_lambda1)
    ADREPORT(sigma_lambda2)
    ADREPORT(sigma_lambdac)

    # annual index (full grid, area-weighted) reuses the same surface
    annual_index <- numeric(n_years)
    for (t in 1:n_years) {
      annual_index[t] <- sum(exp(log_expected_density_pred[, t]) * area_pred)
    }
    log_index <- log(annual_index)
    ADREPORT(annual_index)
    ADREPORT(log_index)

    # variance priors
    log_sigma_O_pos %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_E_pos %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_O_bin %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_E_bin %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_obs %~% dnorm(log(0.5), 0.5, log = TRUE)

    ADREPORT(range_shared)
    ADREPORT(year_effects_pos)
    ADREPORT(year_effects_bin)
  }
}

#############################################################
# Buffer footprint + plug-in index (as in the sim)
#############################################################
create_buffer_grid <- function(sample_data, full_grid, buffer_km = 10) {
  result <- NULL
  for (yr in unique(full_grid$year)) {
    yr_grid <- dplyr::filter(full_grid, year == yr)
    yr_samples <- dplyr::filter(sample_data, year == yr)
    dx <- outer(yr_grid$X, yr_samples$X, "-")
    dy <- outer(yr_grid$Y, yr_samples$Y, "-")
    min_dist <- apply(sqrt(dx^2 + dy^2), 1, min)
    result <- rbind(result, yr_grid[min_dist <= buffer_km, ])
  }
  result
}

predict_annual_index <- function(sdrep, mesh, pred_grid_yr1, years, cell_area) {
  pars <- as.list(sdrep, "Estimate")
  A_pred <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(pred_grid_yr1[, c("X", "Y")]))
  omega_pred_pos <- as.vector(A_pred %*% pars$omega_s_pos)
  omega_pred_bin <- as.vector(A_pred %*% pars$omega_s_bin)
  eps_pred_pos <- A_pred %*% pars$epsilon_st_pos
  eps_pred_bin <- A_pred %*% pars$epsilon_st_bin
  annual_index <- vapply(seq_along(years), function(t) {
    ld <- pars$year_effects_pos[t] + omega_pred_pos + eps_pred_pos[, t]
    le <- pars$year_effects_bin[t] + omega_pred_bin + eps_pred_bin[, t]
    sum(plogis(le) * exp(ld) * cell_area)
  }, numeric(1))
  data.frame(year = years, est = annual_index)
}

#############################################################
# Main loop
#############################################################
random_effects <- c(
  "omega_s_pos", "epsilon_st_pos",
  "omega_s_bin", "epsilon_st_bin", "alpha_catch",
  "omega_eta",
  "omega_lambda1", "omega_lambda2", "omega_lambdac"
)

all_index_results <- list()
all_param_results <- list()
all_alpha_results <- list()

for (scenario_name in scenario_names) {
  dat_clean_list <- readRDS(paste0("simulation_paper/data/", scenario_name, "_dc.rds"))

  for (rep in seq_along(dat_clean_list)) {
    dat_clean <- dat_clean_list[[rep]]

    tmb_data <- make_data(
      dat_clean, mesh, spde, A_pred, area_pred_full,
      nearby_all
    )
    par_list <- make_pars(dat_clean, n_vertices, n_years)

    fit_failed <- FALSE
    tryCatch(
      {
        nll <- make_nll_spatiotemporal(tmb_data)
        obj <- MakeADFun(nll, par_list,
          random = random_effects, silent = TRUE,
          inner.control = list(maxit = 1000, tol10 = 1e-3)
        )
        opt <- nlminb(obj$par, obj$fn, obj$gr,
          control = list(iter.max = 2000, eval.max = 2000, abs.tol = 1e-5)
        )
        for (rerun in 1:2) {
          g <- tryCatch(max(abs(obj$gr(opt$par))), error = function(e) Inf)
          if (opt$convergence == 0 && g < 0.01) break
          opt <- nlminb(opt$par, obj$fn, obj$gr,
            control = list(
              iter.max = 2000, eval.max = 4000,
              rel.tol = 1e-9, abs.tol = 1e-9
            )
          )
        }
        sdrep <- sdreport(obj)
        max_grad <- max(abs(sdrep$gradient.fixed))
        converged <- opt$convergence == 0 && max_grad < 0.01
      },
      error = function(e) fit_failed <<- TRUE
    )
    if (fit_failed) {
      next
    }

    ssum <- summary(sdrep, "report")

    # all parameters here should match names of joint model output
    li <- which(rownames(ssum) == "log_index")
    idx_full <- data.frame(
      year     = years,
      est      = exp(ssum[li, "Estimate"]),
      log_est  = ssum[li, "Estimate"],
      se       = ssum[li, "Std. Error"],
      lwr      = exp(ssum[li, "Estimate"] - 1.96 * ssum[li, "Std. Error"]),
      upr      = exp(ssum[li, "Estimate"] + 1.96 * ssum[li, "Std. Error"]),
      grid     = "full"
    )

    # buffer index: plug-in on the year-1 buffer footprint (NA uncertainty)
    buf_grid <- create_buffer_grid(dat_clean, grid, buffer_km = 10)
    buf_yr1 <- dplyr::filter(buf_grid, year == years[1])
    idx_buf_raw <- predict_annual_index(
      sdrep, mesh, buf_yr1, years,
      buf_yr1$area_km2_WCGBTS
    )
    idx_buf <- data.frame(
      year     = idx_buf_raw$year,
      est      = idx_buf_raw$est,
      log_est  = NA_real_,
      se       = NA_real_,
      lwr      = NA_real_,
      upr      = NA_real_,
      grid     = "buffer"
    )

    for (idx_df in list(idx_full, idx_buf)) {
      idx_df$scenario <- scenario_name
      idx_df$model <- "dc" # so this stacks with naive / joint
      idx_df$rep <- rep
      idx_df$converged <- converged
      all_index_results[[length(all_index_results) + 1L]] <- idx_df
    }

    # DC-specific per-fit parameters in separate table
    rg <- which(rownames(ssum) == "range_shared")
    sa <- which(rownames(ssum) == "sigma_alpha")
    be <- which(rownames(ssum) == "b_eta")
    se <- which(rownames(ssum) == "sigma_eta")
    sl1 <- which(rownames(ssum) == "sigma_lambda1")
    sl2 <- which(rownames(ssum) == "sigma_lambda2")
    slc <- which(rownames(ssum) == "sigma_lambdac")

    all_param_results[[length(all_param_results) + 1L]] <- data.frame(
      scenario = scenario_name,
      model = "dc",
      rep = rep,
      converged = converged,
      range = ssum[rg, "Estimate"],
      sigma_alpha = ssum[sa, "Estimate"],
      b_eta = ssum[be, "Estimate"],
      b_eta_se = ssum[be, "Std. Error"],
      sigma_eta = ssum[se, "Estimate"], sigma_eta_se = ssum[se, "Std. Error"],
      sigma_lambda1 = ssum[sl1, "Estimate"], sigma_lambda1_se = ssum[sl1, "Std. Error"],
      sigma_lambda2 = ssum[sl2, "Estimate"], sigma_lambda2_se = ssum[sl2, "Std. Error"],
      sigma_lambdac = ssum[slc, "Estimate"], sigma_lambdac_se = ssum[slc, "Std. Error"]
    )

    # recovered alpha_catch path (one row per year)
    ap <- which(rownames(ssum) == "alpha_catch")
    all_alpha_results[[length(all_alpha_results) + 1L]] <- data.frame(
      scenario    = scenario_name,
      model       = "dc",
      rep         = rep,
      converged   = converged,
      year        = years,
      alpha_catch = ssum[ap, "Estimate"],
      alpha_se    = ssum[ap, "Std. Error"]
    )
  }
}

#############################################################
# Combine / save: index table matches the joint model; the
# DC-specific parameters and the alpha path go to separate files.
#############################################################
index_results <- dplyr::bind_rows(all_index_results)
param_results <- dplyr::bind_rows(all_param_results)
alpha_results <- dplyr::bind_rows(all_alpha_results)

saveRDS(index_results, paste0("simulation_paper/dchoice_index_", b_tag, "_dc.rds"))
saveRDS(param_results, paste0("simulation_paper/dchoice_params_", b_tag, "_dc.rds"))
saveRDS(alpha_results, paste0("simulation_paper/dchoice_alpha_", b_tag, "_dc.rds"))

rjoint <- readRDS(paste0("simulation_paper/dchoice_params_", "tvb", ".rds"))
rdc <- readRDS(paste0("simulation_paper/dchoice_params_", "tvb", "_dc.rds"))
rjoint$Model <- "Joint"
rdc$Model <- "Discrete-choice"
rall <- rbind(rjoint, rdc)

df1 <- rall[, c("scenario", "model", "rep", "sigma_eta", "Model")] |>
  dplyr::rename(value = "sigma_eta") |>
  dplyr::mutate(par = "sigma_eta")
df2 <- rall[, c("scenario", "model", "rep", "sigma_lambda1", "Model")] |>
  dplyr::rename(value = "sigma_lambda1") |>
  dplyr::mutate(par = "sigma_lambda1")
df3 <- rall[, c("scenario", "model", "rep", "sigma_lambda2", "Model")] |>
  dplyr::rename(value = "sigma_lambda2") |>
  dplyr::mutate(par = "sigma_lambda2")
df4 <- rall[, c("scenario", "model", "rep", "sigma_lambdac", "Model")] |>
  dplyr::rename(value = "sigma_lambdac") |>
  dplyr::mutate(par = "sigma_lambdac")
dfs <- rbind(df1, df2, df3, df4)

ggplot(dfs, aes(value, col = Model, fill = Model)) +
  # geom_vline(aes(xintercept = log(0.5))) +
  geom_histogram(alpha = 0.5, col = NA) +
  facet_grid(par ~ Model) +
  theme_bw() +
  xlab("Estimates from replicates") +
  ylab("Count")
ggsave("simulation_paper/figures/fig_s3_dc_variances.png", width = 7, height = 7)
