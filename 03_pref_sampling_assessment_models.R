library(dplyr)
library(r4ss)

# If true, ss3 models will rerun
# If false, ss3 model results will get read in and written to disk in a nested tibble to help summarize results,
# but models will not run. This is *much* faster.
# Set this to false if you ran the models and encountered issues / made code changes after. Otherwise, true.
run_models <- TRUE

quad_results <- readRDS('simulation_paper/quadrature_results_tvb.rds')

# Create nested data frame, etc.
quad_results_nested <- quad_results |>
  mutate(month = 7, index = 7) |>
  rename(obs = est, se_log = se) |>
  filter(!is.na(se_log)) |>
  select(-lwr, -upr, -log_est) |>
  tidyr::nest(ss3_tab = c(year, obs, month, index, se_log))

# only have windows executable right now, could add a mac-os executable and an OS check if that is important
sst_exe <- here::here('simulation_paper/ss3_models/base_model/ss.exe')

sst_base <- SS_read('simulation_paper/ss3_models/base_model', ss_new = TRUE) # use ss_new to hopefully speed things up

# Add new index fleet. Most of this only needs to be done once, outside of a loop.
sst_add_fleet <- sst_base
# First update data file
sst_add_fleet$dat$Nfleets <- sst_add_fleet$ctl$Nfleets <- sst_base$dat$Nfleets +
  1
sst_add_fleet$dat$fleetinfo <- SS_add_parameter_line(
  sst_base$dat$fleetinfo,
  row_to_copy = 6,
  row_before = 6,
  newval_df = data.frame(fleetname = 'new_index', rowname = 7)
)
sst_add_fleet$dat$fleetnames <- sst_add_fleet$ctl$fleetnames <- sst_add_fleet$dat$fleetinfo$fleetname
sst_add_fleet$dat$CPUEinfo <- SS_add_parameter_line(
  sst_base$dat$CPUEinfo,
  row_to_copy = 6,
  row_before = 6,
  newval_df = data.frame(rowname = 'new_index', fleet = 7)
)
# Even though no new length data, need to add len_info
sst_add_fleet$dat$len_info <- SS_add_parameter_line(
  sst_base$dat$len_info,
  row_to_copy = 6,
  row_before = 6,
  newval_df = data.frame(rowname = 'new_index')
)
# CPUE will get added in the loop, so now move on to control file
# mirror WCGBTS selectivity
sst_add_fleet$ctl$size_selex_types <- bind_rows(
  sst_base$ctl$size_selex_types,
  data.frame(Pattern = c(new_index = 15), Discard = 0, Male = 0, Special = 6) # not sure if male is right
)
# no age selectivity
sst_add_fleet$ctl$age_selex_types <- SS_add_parameter_line(
  sst_base$ctl$age_selex_types,
  row_to_copy = 6,
  row_before = 6,
  newval_df = data.frame(rowname = 'new_index')
)
# Use power function for catchability to see if SS3 alone can account for preferential sampling
# I don't think we should allow down-weighting of the index (i.e., extra_se).
sst_add_fleet$ctl$Q_options <- sst_base$ctl$Q_options |>
  SS_add_parameter_line(
    row_to_copy = 2,
    row_before = 2,
    newval_df = data.frame(
      rowname = 'new_index',
      fleet = 7,
      link = 3,
      float = 0
    )
  )
# Add Q parameters
sst_add_fleet$ctl$Q_parms <- sst_base$ctl$Q_parms |>
  # estimate Q (I believe necessary for a power link function)
  SS_add_parameter_line(
    row_to_copy = 'NWFSCcombo',
    row_before = 'NWFSCcombo',
    newval_df = data.frame(rowname = 'LnQ_new_index', PHASE = 2)
  ) |>
  # start out fixing the power paramter at 0 (i.e., linear relationship)
  SS_add_parameter_line(
    row_to_copy = 'NWFSCcombo',
    row_before = 'new_index',
    newval_df = data.frame(rowname = 'Q_power_new_index', INIT = 0)
  )

future::plan(future::multisession(
  workers = parallelly::availableCores(omit = 1)
))
quad_results_nested_ss3 <- quad_results_nested |>
  mutate(
    ss3_out = furrr::future_imap(ss3_tab, \(df, ii) {
      dir <- file.path('simulation_paper', 'ss3_models', ii)
      if (run_models) {
        sst_temp <- sst_add_fleet
        sst_temp$dat$CPUE <- bind_rows(
          sst_add_fleet$dat$CPUE,
          as.data.frame(df)
        )
        SS_write(sst_temp, dir = dir, overwrite = TRUE)
        r4ss::run(dir, exe = sst_exe, extras = '-nohess', skipfinished = FALSE)
      }
      out <- SS_output(dir = dir, verbose = FALSE, printstats = FALSE)
      return(out[c('timeseries', 'derived_quants', 'parameters', 'maximum_gradient_component')])
    })
  ) # could do this again and allow for estimation of Q_power
future::plan(future::sequential)

unconverged <- quad_results_nested_ss3 |> 
  mutate(max_gradient = purrr::map_dbl(ss3_out, 'maximum_gradient_component')) |> 
  filter(max_gradient > 0.1)

print(paste(nrow(unconverged), 'models did not converge'))

saveRDS(
  quad_results_nested_ss3,
  file = 'simulation_paper/quadrature_ss3_results.rds'
)

# TODO: run SS_plots() on one or two simulations