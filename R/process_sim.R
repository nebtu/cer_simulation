process_fwer <- function(result) {
  result |>
    as_tibble() |>
    group_by(eff, futility, corr, name, run1) |>
    summarise(rej_any = mean(rej_any))
}

process_power <- function(result) {
  num_eff_hyp <- sum(result$eff[[1]] > 0)
  eff_hyp <- paste0(
    "rej_",
    c((5 - num_eff_hyp):4, (9 - num_eff_hyp):8)
  )

  eff_hyp_primary <- paste0(
    "rej_",
    (5 - num_eff_hyp):4
  )

  result |>
    as_tibble() |>
    mutate(
      rej_all_eff = (rowSums(pick(all_of(eff_hyp))) == (num_eff_hyp * 2)),
      rej_all_eff_primary = (rowSums(pick(all_of(eff_hyp_primary))) ==
        num_eff_hyp),
      rej_any_eff = (rowSums(pick(all_of(eff_hyp))) > 0),
      rej_any_eff_primary = (rowSums(pick(all_of(eff_hyp_primary))) > 0)
    ) |>
    group_by(eff, futility, corr, name, run1) |>
    summarise(
      across(starts_with("rej"), mean, .names = "mean_{.col}"),
    )
}

get_power_tbl <- function(power_all_res) {
  tbl_data <- power_all_res |>
    ungroup() |>
    select(
      eff,
      futility,
      mean_rej_all_eff,
      mean_rej_any_eff,
      sd_rej_all_eff,
      sd_rej_any_eff,
      name
    ) |>
    rowwise() |>
    mutate(
      scenario = paste0("S", sum(eff > 0)),
      dropping_rule = case_when(
        futility == 0 ~ "Ultra",
        futility == 0.25 ~ "Aggressive",
        futility == 0.5 ~ "Moderate",
        futility == 0.75 ~ "Conservative"
      ),
      bin = str_ends(name, "bin"),
      conf_low_all = mean_rej_all_eff - (sd_rej_all_eff / sqrt(100000)) * 1.96,
      conf_low_any = mean_rej_any_eff - (sd_rej_any_eff / sqrt(100000)) * 1.96,
      conf_high_all = mean_rej_all_eff + (sd_rej_all_eff / sqrt(100000)) * 1.96,
      conf_high_any = mean_rej_any_eff + (sd_rej_any_eff / sqrt(100000)) * 1.96,
      conf_int_any = paste0(
        "(",
        vec_fmt_percent(conf_low_any, 4),
        ", ",
        vec_fmt_percent(conf_high_any, 4),
        ")"
      ),
      conf_int_all = paste0(
        "(",
        vec_fmt_percent(conf_low_all, 4),
        ", ",
        vec_fmt_percent(conf_high_all, 4),
        ")"
      )
    ) |>
    ungroup() |>
    select(
      scenario,
      dropping_rule,
      mean_rej_any_eff,
      mean_rej_all_eff,
      conf_int_any,
      conf_int_all,
      bin
    ) |>
    group_by(scenario)

  tbl_cont <- tbl_data |>
    filter(!bin) |>
    select(-bin) |>
    gt() |>
    tab_header(
      title = "Disjunctive and Conjunctive Power for continuous scenarios"
    ) |>
    cols_label(
      scenario = "Scenario",
      dropping_rule = "Dropping Rule",
      mean_rej_any_eff = "Disjunctive",
      conf_int_any = "CI (Disj)",
      mean_rej_all_eff = "Conjunctive",
      conf_int_all = "CI (Conj)",
    ) |>
    fmt_percent(
      c(mean_rej_any_eff, mean_rej_all_eff),
      decimals = 4
    )

  tbl_bin <- tbl_data |>
    filter(bin) |>
    select(-bin) |>
    gt() |>
    tab_header(
      title = "Disjunctive and Conjunctive Power for binary scenarios"
    ) |>
    cols_label(
      scenario = "Scenario",
      dropping_rule = "Dropping Rule",
      mean_rej_any_eff = "Disjunctive",
      conf_int_any = "CI (Disj)",
      mean_rej_all_eff = "Conjunctive",
      conf_int_all = "CI (Conj)",
    ) |>
    fmt_percent(
      c(mean_rej_any_eff, mean_rej_all_eff),
      decimals = 4
    )

  list(
    tbl_cont = tbl_cont,
    tbl_bin = tbl_bin
  )
}

get_fwer_tbl <- function(fwer_all_res) {
  tbl_data <- fwer_all_res |>
    mutate(
      dropping_rule = case_when(
        futility == 0 ~ "Ultra",
        futility == 0.25 ~ "Aggressive",
        futility == 0.5 ~ "Moderate",
        futility == 0.75 ~ "Conservative"
      ),
      mean = mean_rej_any,
      bin = str_ends(name, "bin"),
      conf_low = mean - (sd_rej_any / sqrt(100000)) * 1.96,
      conf_high = mean + (sd_rej_any / sqrt(100000)) * 1.96,
      conf_int = paste0(
        "(",
        vec_fmt_percent(conf_low, 4),
        ", ",
        vec_fmt_percent(conf_high, 4),
        ")"
      )
    ) |>
    pivot_wider(
      id_cols = c(corr, bin),
      names_from = dropping_rule,
      values_from = c(mean, conf_int)
    ) |>
    ungroup()

  tbl_bin <- tbl_data |>
    filter(bin) |>
    select(-bin) |>
    gt() |>
    tab_header(
      title = "FWER for binary scenarios"
    ) |>
    cols_label(
      corr = "Correlation",
      mean_Conservative = "Conservative",
      conf_int_Conservative = "CI (Cons)",
      mean_Moderate = "Moderate",
      conf_int_Moderate = "CI (Mod)",
      mean_Aggressive = "Aggressive",
      conf_int_Aggressive = "CI (Aggr)",
      mean_Ultra = "Ultra",
      conf_int_Ultra = "CI (Ultra)"
    ) |>
    fmt_percent(
      mean_Conservative:`conf_int_Ultra`,
      decimals = 4
    )

  tbl_cont <- tbl_data |>
    filter(!bin) |>
    select(-bin) |>
    gt() |>
    tab_header(
      title = "FWER for continuous scenarios"
    ) |>
    cols_label(
      corr = "Correlation",
      mean_Conservative = "Conservative",
      conf_int_Conservative = "CI (Cons)",
      mean_Moderate = "Moderate",
      conf_int_Moderate = "CI (Mod)",
      mean_Aggressive = "Aggressive",
      conf_int_Aggressive = "CI (Aggr)",
      mean_Ultra = "Ultra",
      conf_int_Ultra = "CI (Ultra)"
    ) |>
    fmt_percent(
      mean_Conservative:`conf_int_Ultra`,
      decimals = 4
    )

  list(
    tbl_cont = tbl_cont,
    tbl_bin = tbl_bin
  )
}
