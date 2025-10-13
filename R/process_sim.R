process_fwer <- function(result) {
  result |>
    as_tibble() |>
    group_by(eff, futility, corr, name) |>
    summarise(mean(rej_any), var(rej_any))
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
