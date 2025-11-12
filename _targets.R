# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(future)
library(dplyr)
library(tidyr)
library(autometric)
library(future.callr)

# Set target options:
tar_option_set(
  packages = c(
    "adagraph",
    "future",
    "tidyverse",
    "gt"
  ),
  seed = 1,
  memory = "transient",
  format = "qs",
)
plan(callr)

tar_source(files = "R")

power_analysis <- expand_grid(
  #effect size of primary and secondary arms
  eff = list(
    c(0, 0, 0, 0.4),
    c(0, 0, 0.4, 0.4),
    c(0, 0.4, 0.4, 0.4),
    c(0.4, 0.4, 0.4, 0.4)
  ),
  #correlation between primary and secondary arm
  futility = c(0.75, 0.5, 0.25, 0),
  bin = list(integer(0), c(5, 6, 7, 8))
) |>
  mutate(
    corr = 0.5,
    alt_drop = TRUE,
    runs1 = 1000, #number of completely different trials
    runs2 = 100, #number of second stages simulated per first stage
    n1 = 50 #n2 is automatically the same as n1, since t of the design is 1/2, which is what the data gen function references
  ) |>
  rowwise() |>
  mutate(
    bin_treat_resp = ifelse(
      length(bin) > 0,
      list(ifelse(eff > 0, 0.25, 0.1)),
      list(NULL)
    ),
    bin_con_resp = ifelse(
      length(bin) > 0,
      list(c(0.1, 0.1, 0.1, 0.1)),
      list(NULL)
    )
  ) |>
  mutate(
    eff_name = sum(eff > 0),
    name = paste0(
      "f",
      futility,
      "_e",
      eff_name,
      ifelse(length(bin) > 0, "_bin", "")
    )
  ) |>
  ungroup() |>
  select(-eff_name)

fwer_analysis <- expand_grid(
  corr = c(0, 0.5, 0.8),
  futility = c(0.75, 0.5, 0.25, 0),
  bin = list(integer(0), c(5, 6, 7, 8))
) |>
  mutate(
    eff = list(c(0, 0, 0, 0)),
    alt_drop = TRUE,
    runs1 = 1000, #number of completely different trials
    runs2 = 100, #number of second stages simulated per first stage
    n1 = 50 #n2 is automatically the same as n1, since t of the design is 1/2, which is what the data gen function references
  ) |>
  rowwise() |>
  mutate(
    bin_treat_resp = ifelse(
      length(bin) > 0,
      list(c(0.1, 0.1, 0.1, 0.1)),
      list(NULL)
    ),
    bin_con_resp = ifelse(
      length(bin) > 0,
      list(c(0.1, 0.1, 0.1, 0.1)),
      list(NULL)
    )
  ) |>
  mutate(
    name = paste0(
      "f",
      futility,
      "_c",
      corr,
      ifelse(length(bin) > 0, "_bin", "")
    )
  ) |>
  ungroup()

#power analysis
power_map <- tar_map(
  values = power_analysis,
  names = all_of("name"),
  tar_target(
    index_batch,
    seq_len(100)
  ),
  tar_target(
    data_gen,
    get_sim_data_gen(
      corr = corr,
      eff = eff,
      n1 = n1,
      bin = bin,
      bin_con_resp = bin_con_resp,
      bin_treat_resp = bin_treat_resp
    ),
    deployment = "main"
  ),
  tar_target(
    adaption_func,
    get_sim_adaption(futility = futility, alt_drop = alt_drop),
    deployment = "main"
  ),
  tar_target(
    power,
    sim_trial(
      design,
      runs1 = runs1,
      runs2 = runs2,
      adapt_rule = adaption_func,
      data_gen_1 = data_gen[[1]],
      data_gen_2 = data_gen[[2]]
    ),
    pattern = map(index_batch)
  ),
  tar_target(
    power_summary_batchwise,
    power |>
      as_tibble() |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility,
        name = name
      ) |>
      process_power(),
    pattern = map(power)
  ),
  tar_target(
    power_summary,
    power_summary_batchwise |>
      group_by(across(all_of(c("name", "eff", "corr", "futility")))) |>
      summarise(
        across(
          where(is.numeric),
          list(mean = mean, sd = sd),
          .names = "{.fn}_{.col}"
        ),
      ) |>
      rename_with(
        \(x) str_replace(x, "mean_mean", "mean"),
        .cols = starts_with("mean_mean")
      ) |>
      rename_with(
        \(x) str_replace(x, "sd_mean", "sd"),
        .cols = starts_with("sd_mean")
      ),
    deployment = "main"
  )
)

#fwer analysis
fwer_map <- tar_map(
  values = fwer_analysis,
  names = all_of("name"),
  tar_target(
    fwer_index_batch,
    seq_len(100)
  ),
  tar_target(
    fwer_data_gen,
    get_sim_data_gen(
      corr = corr,
      eff = eff,
      n1 = n1,
      bin = bin,
      bin_con_resp = bin_con_resp,
      bin_treat_resp = bin_treat_resp
    ),
    deployment = "main"
  ),
  tar_target(
    fwer_adaption_func,
    get_sim_adaption(futility = futility, alt_drop = alt_drop),
    deployment = "main"
  ),
  tar_target(
    fwer,
    sim_trial(
      design,
      runs1 = runs1,
      runs2 = runs2,
      adapt_rule = fwer_adaption_func,
      data_gen_1 = fwer_data_gen[[1]],
      data_gen_2 = fwer_data_gen[[2]]
    ),
    pattern = map(fwer_index_batch)
  ),
  tar_target(
    fwer_summary_batchwise,
    fwer |>
      as_tibble() |>
      mutate(
        corr = corr,
        eff = list(eff),
        futility = futility,
        name = name
      ) |>
      process_fwer(),
    pattern = map(fwer)
  ),
  tar_target(
    fwer_summary,
    fwer_summary_batchwise |>
      group_by(across(all_of(c("name", "eff", "corr", "futility")))) |>
      summarise(
        across(
          where(is.numeric),
          list(mean = mean, sd = sd),
          .names = "{.fn}_{.col}"
        ),
      ),
    deployment = "main"
  )
)

list(
  tar_target(
    design,
    get_sim_design(100, 0.5)
  ),
  power_map,
  fwer_map,
  tar_combine(power_all_res, power_map["power_summary"]),
  tar_combine(fwer_all_res, fwer_map["fwer_summary"]),
  tar_target(
    power_tbl,
    get_power_tbl(power_all_res)
  ),
  tar_target(
    fwer_tbl,
    get_fwer_tbl(fwer_all_res)
  ),
  tar_quarto(
    binary_report,
    "binary.qmd"
  )
)
