process_fwer <- function(result) {
  result |>
    as_tibble() |>
    group_by(eff, futility, corr, name) |>
    summarise(mean(rej_any))
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
    group_by(eff, futility, corr, name) |>
    summarise(
      across(starts_with("rej"), mean, .names = "mean_{.col}"),
      across(
        starts_with(c("rej_all", "rej_any")),
        \(x) confint(lm(x ~ 1)),
        .names = "conf_{.col}"
      )
    )
}

get_power_tbl <- function(power_all_res) {
  tbl_data <- power_all_res |>
    ungroup() |>
    select(eff, futility, mean_rej_all_eff, mean_rej_any_eff, name) |>
    rowwise() |>
    mutate(
      scenario = paste0("S", sum(eff > 0)),
      dropping_rule = case_when(
        futility == 0 ~ "Ultra Aggressive",
        futility == 0.25 ~ "Aggressive",
        futility == 0.5 ~ "Moderate",
        futility == 0.75 ~ "Conservative"
      ),
      bin = str_ends(name, "bin")
    ) |>
    ungroup() |>
    select(scenario, dropping_rule, mean_rej_any_eff, mean_rej_all_eff, bin) |>
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
      mean_rej_all_eff = "Conjunctive"
    ) |>
    fmt_number(
      c(mean_rej_any_eff, mean_rej_all_eff),
      decimals = 3
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
      mean_rej_all_eff = "Conjunctive"
    ) |>
    fmt_number(
      c(mean_rej_any_eff, mean_rej_all_eff),
      decimals = 3
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
        futility == 0 ~ "Ultra Aggressive",
        futility == 0.25 ~ "Aggressive",
        futility == 0.5 ~ "Moderate",
        futility == 0.75 ~ "Conservative"
      ),
      mean = `mean(rej_any)`,
      bin = str_ends(name, "bin")
    ) |>
    pivot_wider(
      id_cols = c(corr, bin),
      names_from = dropping_rule,
      values_from = mean
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
      corr = "Correlation"
    ) |>
    fmt_number(
      Conservative:`Ultra Aggressive`,
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
      corr = "Correlation"
    ) |>
    fmt_number(
      Conservative:`Ultra Aggressive`,
      decimals = 4
    )

  list(
    tbl_cont = tbl_cont,
    tbl_bin = tbl_bin
  )
}
