library(dplyr)
library(r4ss)


# only have windows executable right now, could add a mac-os executable and an OS check if that is important
sst_exe <- here::here('simulation_paper/ss3_models/base_model/ss.exe')

sst_base <- SS_read('simulation_paper/ss3_models/base_model')
# This model is slow. I wonder if it is the discards and if those are worth simplifying.
# test if I can add discards to the catch instead of estimating retention:
out <- SS_output('simulation_paper/ss3_models/base_model')
sst_discards <- sst_base
# Remove discard data
sst_discards$dat$discard_fleet_info <- NULL
sst_discards$dat$discard_data <- NULL
sst_discards$dat$N_discard_fleets <- 0
sst_discards$dat$meanbodywt <- NULL
sst_discards$dat$use_meanbodywt <- 0
sst_discards$dat$DF_for_meanbodywt <- NULL
sst_discards$dat$lencomp <- sst_base$dat$lencomp |>
  filter(part != 1) |> # exclude discard comps
  mutate(part = 0) # make remaining comps combined discards + retained
# this will set catch = total estimated dead biomass (by fleet) from the model w/retention
sst_discards$dat$catch <- sst_base$dat$catch |>
  left_join(
    select(out$catch, Fleet, Yr, kill_bio),
    by = c(year = 'Yr', fleet = 'Fleet')
  ) |>
  select(-catch) |>
  rename(catch = kill_bio) |>
  filter(!is.na(catch)) # there is an NA in 2023 this excludes
# remove retention parameters
sst_discards$ctl$size_selex_types$Discard <- 0
sst_discards$ctl$size_selex_parms <- sst_base$ctl$size_selex_parms[
  -grep('Ret', rownames(sst_base$ctl$size_selex_parms)),
]
sst_discards$ctl$size_selex_parms_tv <- sst_base$ctl$size_selex_parms_tv[
  -grep('Ret', rownames(sst_base$ctl$size_selex_parms_tv)),
]
SS_write(
  sst_discards,
  'simulation_paper/ss3_models/no_retention',
  overwrite = TRUE
)
run(
  'simulation_paper/ss3_models/no_retention',
  exe = sst_exe,
  extras = '-nohess',
  show_in_console = TRUE
)
# this was indeed MUCH faster (10:12 vs 1:17).
# now see if things are appreciably different
out_list <- r4ss::SSgetoutput(
  dirvec = c(
    'simulation_paper/ss3_models/no_retention',
    'simulation_paper/ss3_models/base_model'
  ),
  modelnames = c('no_retention', 'base_2023')
)
out_list |>
  SSsummarize() |>
  SSplotComparisons(subplots = c(1, 3, 11), new = FALSE)
# They are extremely different. :(
# Especially recdevs before ~2000, which propagates to unfished stuff
# Is this a situation like widow where they are different but shouldn't be?
# Discarding was commonplace in those years so maybe it's ok they are different.
# Looking at fits to length comps, it looks like for the trawl fleet in particular
# the discarded fish are indeed much smaller than the retained fish.
# It might be reasonable to at least remove retention on the non-trawl fleet.
# One question is: can I still use this (overly) simplified model for this purpose?
# Another alternative is to set up discard fleets. That would probably also run faster, but more work.

# Let's instead fix the retention parameters
# rerun model to get .ss_new files
run(
  'simulation_paper/ss3_models/base_model',
  exe = sst_exe,
  extras = '-nohess',
  show_in_console = FALSE,
  skipfinished = FALSE
) # actually on my new laptop this only took 4:12
mod_ss_new <- SS_read('simulation_paper/ss3_models/base_model', ss_new = TRUE)
mod_ss_new$ctl$size_selex_parms[
  grep('Ret', rownames(mod_ss_new$ctl$size_selex_parms)),
  'PHASE'
] <- -99
mod_ss_new$ctl$size_selex_parms_tv[
  grep('Ret', rownames(mod_ss_new$ctl$size_selex_parms_tv)),
  'PHASE'
] <- -99
mod_ss_new$dat$lencomp <- mod_ss_new$dat$lencomp |>
  filter(part != 1)
mod_ss_new$dat$use_meanbodywt <- 0
mod_ss_new$dat$DF_for_meanbodywt <- NULL
mod_ss_new$dat$meanbodywt <- NULL
SS_write(
  mod_ss_new,
  'simulation_paper/ss3_models/fix_retention',
  overwrite = TRUE
)
run(
  'simulation_paper/ss3_models/fix_retention',
  exe = sst_exe,
  extras = '-nohess',
  skipfinished = FALSE
) # 3:33. I give up. I am leaving it as is.
