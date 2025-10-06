# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(future)
library(dplyr)
library(tidyr)


# Set target options:
tar_option_set(
  packages = c(
    "adagraph",
    "future",
    "tidyverse"
  ),
  seed = 1,
  memory = "transient",
  garbage_collection = 2
)

tar_source(files = "R")
plan(multicore)

power_analysis <- expand_grid(
  #effect size of primary and secondary arms
  eff = list(
    c(0, 0, 0, 0.4),
    c(0, 0, 0.4, 0.4),
    c(0, 0.4, 0.4, 0.4),
    c(0.4, 0.4, 0.4, 0.4)
  ),
  #correlation between primary and secondary arm
  futility = c(0.75, 0.5, 0.25, 0)
) |>
  mutate(corr = 0.5) |>
  rowwise() |>
  mutate(eff_name = sum(eff > 0), name = paste0(futility, "_", eff_name)) |>
  ungroup() |>
  select(-eff_name)

power_analysis_bin <- power_analysis |>
  rowwise() |>
  mutate(
    bin_prop = list(ifelse(eff > 0, 0.25, 0.1))
  ) |>
  ungroup()

fwer_analysis <- expand_grid(
  corr = c(0, 0.5, 0.8),
  futility = c(0.75, 0.5, 0.25, 0)
) |>
  mutate(eff = list(c(0, 0, 0, 0))) |>
  mutate(name = paste0("f", futility, "_", "c", corr))


#power analysis
power_map <- tar_map(
  values = power_analysis,
  names = any_of("name"),
  tar_target(
    power,
    run_example_trial(
      design,
      runs1 = 1000,
      runs2 = 100,
      n1 = 50,
      n2 = 50,
      cor = corr,
      eff = eff,
      futility = futility
    )
  ),
  tar_target(
    power_summary,
    power |>
      as_tibble() |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility
      ) |>
      process_power()
  )
)

power_map_bin <- tar_map(
  values = power_analysis_bin,
  names = any_of("name"),
  tar_target(
    power_bin,
    run_example_trial_bin(
      design,
      runs1 = 1000,
      runs2 = 100,
      n1 = 50,
      n2 = 50,
      cor = corr,
      eff = eff,
      futility = futility,
      bin_con_resp = c(0.1, 0.1, 0.1, 0.1),
      bin_treat_resp = bin_prop
    )
  ),
  tar_target(
    power_summary_bin,
    power_bin |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility
      ) |>
      process_power()
  )
)

#fwer_analysis
fwer_map <- tar_map(
  values = fwer_analysis,
  names = any_of("name"),
  tar_target(
    fwer,
    run_example_trial(
      design,
      runs1 = 1000,
      runs2 = 100,
      n1 = 50,
      n2 = 50,
      cor = corr,
      eff = eff,
      futility = futility
    )
  ),
  tar_target(
    fwer_summary,
    fwer |>
      as_tibble() |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility
      ) |>
      process_fwer()
  )
)

fwer_map_bin <- tar_map(
  values = fwer_analysis,
  names = any_of("name"),
  tar_target(
    fwer_bin,
    run_example_trial_bin(
      design,
      runs1 = 1000,
      runs2 = 100,
      n1 = 50,
      n2 = 50,
      cor = corr,
      eff = eff,
      futility = futility,
      bin_con_resp = c(0.1, 0.1, 0.1, 0.1),
      bin_treat_resp = c(0.1, 0.1, 0.1, 0.1)
    )
  ),
  tar_target(
    fwer_summary_bin,
    fwer_bin |>
      as_tibble() |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility
      ) |>
      process_fwer()
  )
)
list(
  tar_target(
    design,
    get_sim_design(100, 0.5)
  ),
  power_map,
  fwer_map,
  power_map_bin,
  fwer_map_bin,
  tar_combine(power_all_res, power_map["power_summary"]),
  tar_combine(fwer_all_res, fwer_map["fwer_summary"]),
  tar_combine(power_all_res_bin, power_map_bin["power_summary_bin"]),
  tar_combine(fwer_all_res_bin, fwer_map_bin["fwer_summary_bin"])
)
