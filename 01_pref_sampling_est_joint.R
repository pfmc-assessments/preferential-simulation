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

# Time-varying preferential sampling:
#   If time_varying_b = TRUE, the preferential coefficient follows a linear
#   trend (0 in year 1 and increases to b_max in the final year) plus AR(1) process
#   deviations. If FALSE, b is CONSTANT at b_max in every year (level-bias
#   regime). b_max indexes the scenario; b_max = 0 is the random control.
time_varying_b <- TRUE # FALSE = constant b
dc_model <- FALSE
b_max_grid <- c(0, 3, 5) # c(0, 1, 2, 3, 5)

# scenario list
scenarios <- setNames(as.list(b_max_grid), paste0("b", b_max_grid))

#############################################################
# Buffered grid (prediction footprint shrinks under strong pref sampling)
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

#############################################################
# Build the sampling indicator matrix R (n_grid x n_years)
#############################################################
make_R_matrix <- function(dat_clean, n_grid, years) {
  R <- matrix(0L, n_grid, length(years))
  yi <- match(dat_clean$year, years)
  R[cbind(dat_clean$cell_id, yi)] <- 1L
  R
}

#############################################################
# RTMB data list
#   use_rw_b: 0 = scalar b (joint), 1 = random-walk b_t (joint_rw)
#############################################################
make_data <- function(dat_clean, mesh, spde, A_grid, R_mat, use_preferential,
                      area_pred, use_rw_b = 0L) {
  years_v <- sort(unique(dat_clean$year))
  year_idx <- match(dat_clean$year, years_v)
  A_obs <- fmesher::fm_basis(
    mesh$mesh,
    loc = as.matrix(dat_clean[, c("X", "Y")])
  )
  list(
    observations = dat_clean$total_catch_wt_kg,
    log_effort = dat_clean$log_effort,
    year_indices = year_idx,
    pos_idx = which(dat_clean$total_catch_wt_kg > 0),
    n_obs = nrow(dat_clean),
    n_years = length(years_v),
    use_preferential = as.numeric(use_preferential),
    use_rw_b = as.numeric(use_rw_b), # 0 = scalar b, 1 = random walk b_t
    spde = spde,
    A_obs = A_obs,
    A_grid = A_grid, # mesh -> all grid cells (also used as A_pred)
    R_mat = R_mat, # sampling indicator (n_grid x n_years)
    n_grid = nrow(A_grid),
    n_vertices = mesh$mesh$n,
    area_pred = area_pred # full-grid cell areas for the index
  )
}


#############################################################
# RTMB parameter list
#############################################################
make_pars <- function(dat_clean, n_vertices, n_years, use_preferential,
                      use_rw_b = 0L) {
  par <- list(
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
    epsilon_st_bin = matrix(0, n_vertices, n_years)
  )
  if (use_preferential == 1L) {
    par$gamma_0 <- 0
    par$eta_s <- rep(0, n_vertices)
    par$log_sigma_eta <- log(0.5)
    if (use_rw_b == 1L) {
      par$b_pref <- rep(0, n_years) # Conn b as a random-walk path
      par$log_sigma_b <- log(0.3) # RW innovation SD hyperparameter
    } else {
      par$b_pref <- 0 # Conn b (scalar, unconstrained)
    }
  }
  par
}


#############################################################
# NLL - main function
#############################################################
make_joint_nll_spatiotemporal <- function(tmb_data) {
  function(par) {
    getAll(par, tmb_data)
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    kappa <- exp(log_kappa)
    range_shared <- sqrt(8) / kappa
    sigma_O_pos <- exp(log_sigma_O_pos)
    sigma_E_pos <- exp(log_sigma_E_pos)
    sigma_O_bin <- exp(log_sigma_O_bin)
    sigma_E_bin <- exp(log_sigma_E_bin)
    sigma_obs <- exp(log_sigma_obs)

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

    # process model (delta-lognormal) at the hauls
    omega_proj_pos <- as.vector(A_obs %*% omega_s_pos)
    omega_proj_bin <- as.vector(A_obs %*% omega_s_bin)
    eps_proj_pos <- (A_obs %*% epsilon_st_pos)[cbind(seq_len(n_obs), year_indices)]
    eps_proj_bin <- (A_obs %*% epsilon_st_bin)[cbind(seq_len(n_obs), year_indices)]

    log_density_pos <- year_effects_pos[year_indices] +
      omega_proj_pos + eps_proj_pos
    logit_encounter <- year_effects_bin[year_indices] +
      omega_proj_bin + eps_proj_bin
    log_abundance_pos <- log_density_pos + log_effort

    as.numeric(observations > 0) %~%
      dbinom(size = 1, prob = plogis(logit_encounter))
    log(observations[pos_idx]) %~%
      dnorm(log_abundance_pos[pos_idx], sigma_obs)

    # Conn-style areal sampling model
    if (use_preferential == 1) {
      # independent sampling field (eta)
      sigma_eta <- exp(log_sigma_eta)
      tau_eta <- 1 / (sqrt(4 * pi) * kappa * sigma_eta)
      Q_eta <- Q_base * tau_eta^2
      eta_s %~% dgmrf(0, Q_eta)
      eta_grid <- as.vector(A_grid %*% eta_s)

      # abundance spatial field deviation at every grid cell (Conn delta)
      omega_grid_pos <- as.vector(A_grid %*% omega_s_pos)
      eps_grid_pos <- A_grid %*% epsilon_st_pos # n_grid x n_years

      if (use_rw_b == 1) {
        # time-varying b_t via random walk
        sigma_b_rw <- exp(log_sigma_b)
        # diffuse prior on year-1; RW for subsequent years
        b_pref[1] %~% dnorm(0, 2.0, log = TRUE)
        for (t in 2:n_years) {
          b_pref[t] %~% dnorm(b_pref[t - 1], sigma_b_rw, log = TRUE)
        }
        # Bernoulli over the full grid each year with the years own b_t
        for (t in 1:n_years) {
          delta_t <- omega_grid_pos + eps_grid_pos[, t]
          logit_R_t <- gamma_0 + b_pref[t] * delta_t + eta_grid
          R_mat[, t] %~% dbinom(size = 1, prob = plogis(logit_R_t))
        }
        log_sigma_b %~% dnorm(log(0.3), 0.5, log = TRUE) # RW innovation SD
        ADREPORT(sigma_b_rw)
        ADREPORT(b_pref) # length-n_years path (with SEs)
      } else {
        # scalar b (Conn model)
        for (t in 1:n_years) {
          delta_t <- omega_grid_pos + eps_grid_pos[, t]
          logit_R_t <- gamma_0 + b_pref * delta_t + eta_grid
          R_mat[, t] %~% dbinom(size = 1, prob = plogis(logit_R_t))
        }
        b_pref %~% dnorm(0, 2.0, log = TRUE)
        ADREPORT(b_pref)
      }

      gamma_0 %~% dnorm(0, 1.0, log = TRUE)
      log_sigma_eta %~% dnorm(log(0.5), 0.5, log = TRUE)
      ADREPORT(sigma_eta)
      ADREPORT(gamma_0)
    }

    # Annual index as a derived quantity
    # Full-grid prediction (A_grid projects mesh to grid_year1 cells).
    # Reported on the log scale: SE ~ CV
    omega_pred_pos <- as.vector(A_grid %*% omega_s_pos)
    omega_pred_bin <- as.vector(A_grid %*% omega_s_bin)
    eps_pred_pos <- A_grid %*% epsilon_st_pos
    eps_pred_bin <- A_grid %*% epsilon_st_bin

    log_index <- numeric(n_years)
    for (t in 1:n_years) {
      ld <- year_effects_pos[t] + omega_pred_pos + eps_pred_pos[, t]
      le <- year_effects_bin[t] + omega_pred_bin + eps_pred_bin[, t]
      log_index[t] <- log(sum(plogis(le) * exp(ld) * area_pred))
    }
    ADREPORT(log_index)

    # variance priors
    log_sigma_O_pos %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_E_pos %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_O_bin %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_E_bin %~% dnorm(log(0.5), 0.5, log = TRUE)
    log_sigma_obs %~% dnorm(log(0.5), 0.5, log = TRUE)

    ADREPORT(range_shared)
    ADREPORT(sigma_O_pos)
    ADREPORT(sigma_E_pos)
    ADREPORT(sigma_O_bin)
    ADREPORT(sigma_E_bin)
    ADREPORT(sigma_obs)
    ADREPORT(year_effects_pos)
    ADREPORT(year_effects_bin)
  }
}


#############################################################
# Index prediction (uses the prediction grid's own cell areas)
#############################################################
predict_annual_index <- function(sdrep, mesh, pred_grid_yr1, years, cell_area) {
  pars <- as.list(sdrep, "Estimate")
  A_pred <- fmesher::fm_basis(
    mesh$mesh,
    loc = as.matrix(pred_grid_yr1[, c("X", "Y")])
  )
  omega_pred_pos <- as.vector(A_pred %*% pars$omega_s_pos)
  omega_pred_bin <- as.vector(A_pred %*% pars$omega_s_bin)
  eps_pred_pos <- A_pred %*% pars$epsilon_st_pos
  eps_pred_bin <- A_pred %*% pars$epsilon_st_bin

  annual_index <- vapply(seq_along(years), function(t) {
    log_density_pred <- pars$year_effects_pos[t] + omega_pred_pos + eps_pred_pos[, t]
    logit_enc_pred <- pars$year_effects_bin[t] + omega_pred_bin + eps_pred_bin[, t]
    sum(plogis(logit_enc_pred) * exp(log_density_pred) * cell_area)
  }, numeric(1))

  data.frame(year = years, est = annual_index)
}


#############################################################
# Main loop
#############################################################

random_effects <- c(
  "omega_s_pos", "epsilon_st_pos",
  "omega_s_bin", "epsilon_st_bin"
)

all_index_results <- list()
all_param_results <- list()
all_btraj_results <- list()
all_bpath_results <- list() # RW b_t path per (scenario, rep, year)
result_counter <- 1L

# model spec table: label (use_preferential, use_rw_b)
#   naive : no sampling model
#   joint : sampling model with random-walk b_t
model_table <- list(
  naive = c(use_pref = 0L, use_rw = 0L),
  joint = c(use_pref = 1L, use_rw = 1L) # Conn areal model with RW b_t
)

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

# estimation mesh (cutoff 20km) + derived objects
mesh <- sdmTMB::make_mesh(grid_year1, xy_cols = c("X", "Y"), cutoff = 20)
spde <- fmesher::fm_fem(mesh$mesh)
n_vertices <- mesh$mesh$n
A_grid <- fmesher::fm_basis(mesh$mesh, loc = as.matrix(grid_year1[, c("X", "Y")]))
area_pred_full <- grid_year1$area_km2_WCGBTS


for (scenario_name in names(scenarios)) {
  if (dc_model == TRUE) {
    dat_clean_list <- readRDS(file = paste0("simulation_paper/data/", scenario_name, "_dc.rds"))
  } else {
    dat_clean_list <- readRDS(file = paste0("simulation_paper/data/", scenario_name, ".rds"))
  }

  n_reps <- length(dat_clean_list)
  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100)

    dat_clean <- dat_clean_list[[rep]]
    R_mat <- make_R_matrix(dat_clean, n_grid, years)

    # save this replicate's realized b_t trajectory (true time-varying b)
    all_btraj_results[[length(all_btraj_results) + 1L]] <- data.frame(
      scenario = scenario_name, rep = rep, year = years,
      b_true = as.numeric(attr(dat_clean, "b_t"))
    )

    buf_grid <- create_buffer_grid(dat_clean, grid, buffer_km = 10)
    buf_yr1 <- dplyr::filter(buf_grid, year == years[1])

    # b0 is the random-sampling control: only the naive model is meaningful.
    model_labels <- if (scenario_name == "b0") "naive" else names(model_table)

    for (model_label in model_labels) {
      use_pref <- model_table[[model_label]]["use_pref"]
      use_rw <- model_table[[model_label]]["use_rw"]

      tmb_data <- make_data(
        dat_clean, mesh, spde, A_grid, R_mat, use_pref,
        area_pred_full,
        use_rw_b = use_rw
      )
      par_list <- make_pars(dat_clean, n_vertices, n_years, use_pref,
        use_rw_b = use_rw
      )

      fit_failed <- FALSE
      tryCatch(
        {
          nll <- make_joint_nll_spatiotemporal(tmb_data)
          re <- random_effects
          if (use_pref == 1L) re <- c(re, "eta_s")
          if (use_rw == 1L) re <- c(re, "b_pref")

          obj <- MakeADFun(
            nll, par_list,
            random = re, silent = TRUE,
            inner.control = list(maxit = 1000, tol10 = 1e-3)
          )

          opt <- nlminb(
            obj$par, obj$fn, obj$gr,
            control = list(iter.max = 2000, eval.max = 2000, abs.tol = 1e-5)
          )
          for (rerun in 1:2) {
            g <- tryCatch(max(abs(obj$gr(opt$par))), error = function(e) Inf)
            if (opt$convergence == 0 && g < 0.01) break
            opt <- nlminb(
              opt$par, obj$fn, obj$gr,
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
      if (fit_failed) next

      # Annual index for each grid
      # full grid: read the ADREPORT log_index (with SE);
      # buffer grid: plug-in point estimate
      ssum <- summary(sdrep, "report")
      li_rows <- which(rownames(ssum) == "log_index")

      idx_full <- data.frame(
        year    = years,
        est     = exp(ssum[li_rows, "Estimate"]),
        log_est = ssum[li_rows, "Estimate"],
        se      = ssum[li_rows, "Std. Error"],
        lwr     = exp(ssum[li_rows, "Estimate"] - 1.96 * ssum[li_rows, "Std. Error"]),
        upr     = exp(ssum[li_rows, "Estimate"] + 1.96 * ssum[li_rows, "Std. Error"]),
        grid    = "full"
      )

      idx_buf <- predict_annual_index(
        sdrep, mesh, buf_yr1, years,
        buf_yr1$area_km2_WCGBTS
      )
      idx_buf$log_est <- NA_real_
      idx_buf$se <- NA_real_
      idx_buf$lwr <- NA_real_
      idx_buf$upr <- NA_real_
      idx_buf$grid <- "buffer"

      for (idx_df in list(idx_full, idx_buf)) {
        idx_df$scenario <- scenario_name
        idx_df$model <- model_label
        idx_df$rep <- rep
        idx_df$converged <- converged
        all_index_results[[result_counter]] <- idx_df
        result_counter <- result_counter + 1L
      }

      # get parameters
      p_fixed <- sdrep$par.fixed
      g0 <- if ("gamma_0" %in% names(p_fixed)) p_fixed["gamma_0"] else NA_real_
      lk <- if ("log_kappa" %in% names(p_fixed)) p_fixed["log_kappa"] else NA_real_

      bp_scalar <- NA_real_ # b summary (scalar model = b; RW model = path mean)
      sigma_b_val <- NA_real_ # RW innovation SD (RW joint only)

      if (use_rw == 1L) {
        # time-varying b: RW path is ADREPORT (n_years rows named "b_pref")
        bp_rows <- which(rownames(ssum) == "b_pref")
        if (length(bp_rows) == n_years) {
          all_bpath_results[[length(all_bpath_results) + 1L]] <- data.frame(
            scenario  = scenario_name,
            model     = model_label,
            rep       = rep,
            converged = converged,
            year      = years,
            b_pref    = ssum[bp_rows, "Estimate"],
            b_pref_se = ssum[bp_rows, "Std. Error"]
          )
          bp_scalar <- mean(ssum[bp_rows, "Estimate"])
        }
        sb_rows <- which(rownames(ssum) == "sigma_b_rw")
        if (length(sb_rows) == 1L) sigma_b_val <- ssum[sb_rows, "Estimate"]
      } else if (use_pref == 1L) {
        # scalar b is a fixed effect
        bp_scalar <- if ("b_pref" %in% names(p_fixed)) {
          as.numeric(p_fixed["b_pref"])
        } else {
          NA_real_
        }
      }

      all_param_results[[length(all_param_results) + 1L]] <- data.frame(
        scenario  = scenario_name,
        model     = model_label,
        rep       = rep,
        converged = converged,
        range     = if (!is.na(lk)) sqrt(8) / exp(as.numeric(lk)) else NA_real_,
        gamma_0   = as.numeric(g0),
        b_pref    = bp_scalar, # scalar b, or RW path mean
        sigma_b   = sigma_b_val # RW innovation SD (NA for naive/joint)
      )
    }
  }
}


#############################################################
# Combine / save
#############################################################

# regime for output filenames
b_tag <- if (time_varying_b) "tvb" else "constb"
dc_tag <- if (dc_model) "_dc" else ""
results <- dplyr::bind_rows(all_index_results)
saveRDS(results, paste0("simulation_paper/quadrature_results_", b_tag, dc_tag, ".rds"))

results <- results |>
  dplyr::group_by(scenario, model, grid, rep) |>
  dplyr::mutate(scaled_est = est / est[year == min(year)]) |>
  dplyr::ungroup()

truth_df <- data.frame(
  year = index_baseline$year, est = index_baseline$est,
  scenario = "truth", model = "sdmTMB", grid = "full", rep = 0, converged = TRUE
) |> dplyr::mutate(scaled_est = est / est[year == min(year)])

results_all <- dplyr::bind_rows(results, truth_df)
params <- dplyr::bind_rows(all_param_results) |> dplyr::filter(converged)
b_traj_df <- dplyr::bind_rows(all_btraj_results)
b_pref_path <- dplyr::bind_rows(all_bpath_results) # RW estimated b_t paths

b0_ref <- results_all |>
  dplyr::filter(scenario == "b0", model == "naive") |>
  dplyr::select(year, grid, rep, b0_est = est)
results_all <- results_all |>
  dplyr::left_join(b0_ref, by = c("year", "grid", "rep")) |>
  dplyr::mutate(index_rel_b0 = est / b0_est)

saveRDS(
  list(
    results = results_all, params = params,
    truth = truth_df, b_traj = b_traj_df,
    b_pref_path = b_pref_path
  ),
  paste0("simulation_paper/pref_sampling_results_", b_tag, dc_tag, ".rds")
)
