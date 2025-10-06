#'@importFrom stats pnorm qnorm pt
library(adagraph)

get_sim_design <- function(n, t) {
  ws <- .75 # (weight given to secondary endpoint)
  wp <- (1 - ws) / 3 # weight given to the other arms primary endpoint
  m <- rbind(
    H1 = c(0, wp, wp, wp, ws, 0, 0, 0),
    H2 = c(wp, 0, wp, wp, 0, ws, 0, 0),
    H3 = c(wp, wp, 0, wp, 0, 0, ws, 0),
    H4 = c(wp, wp, wp, 0, 0, 0, 0, ws),
    H5 = c(0, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0),
    H6 = c(1 / 3, 0, 1 / 3, 1 / 3, 0, 0, 0, 0),
    H7 = c(1 / 3, 1 / 3, 0, 1 / 3, 0, 0, 0, 0),
    H8 = c(1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0)
  )

  weights <- c(rep(1 / 4, 4), rep(0, 4))

  alpha <- 0.025 #overall alpha

  #spending function
  as <- function(x, t) 2 - 2 * pnorm(qnorm(1 - x / 2) / sqrt(t))

  design <- multiarm_cer_design(
    controls = 2,
    treatment_assoc = c(1, 1, 1, 1, 2, 2, 2, 2),
    n_controls = n,
    n_treatments = n,
    weights = weights,
    t = t,
    alpha = alpha,
    test_m = m,
    alpha_spending_f = as
  )

  design
}


example_data_gen <- function(corr, eff, n1) {
  mu <- c(eff, eff)

  corr_control <- rbind(
    c(1, corr),
    c(corr, 1)
  )

  corr_treatment <- rbind(
    c(1, 0, 0, 0, corr, 0, 0, 0),
    c(0, 1, 0, 0, 0, corr, 0, 0),
    c(0, 0, 1, 0, 0, 0, corr, 0),
    c(0, 0, 0, 1, 0, 0, 0, corr),
    c(corr, 0, 0, 0, 1, 0, 0, 0),
    c(0, corr, 0, 0, 0, 1, 0, 0),
    c(0, 0, corr, 0, 0, 0, 1, 0),
    c(0, 0, 0, corr, 0, 0, 0, 1)
  )

  data_gen_1 <- get_data_gen(
    corr_control,
    corr_treatment,
    mu,
    n1,
    n1
  )
  data_gen_2 <- get_data_gen_2(
    corr_control,
    corr_treatment,
    mu
  )

  c(data_gen_1, data_gen_2)
}

example_data_gen_bin <- function(
  corr,
  eff,
  n1,
  bin_con_resp,
  bin_treat_resp
) {
  mu <- c(eff, c(NA, NA, NA, NA))

  corr_control <- rbind(
    c(1, corr),
    c(corr, 1)
  )

  corr_treatment <- rbind(
    c(1, 0, 0, 0, corr, 0, 0, 0),
    c(0, 1, 0, 0, 0, corr, 0, 0),
    c(0, 0, 1, 0, 0, 0, corr, 0),
    c(0, 0, 0, 1, 0, 0, 0, corr),
    c(corr, 0, 0, 0, 1, 0, 0, 0),
    c(0, corr, 0, 0, 0, 1, 0, 0),
    c(0, 0, corr, 0, 0, 0, 1, 0),
    c(0, 0, 0, corr, 0, 0, 0, 1)
  )

  data_gen_1 <- get_data_gen(
    corr_control,
    corr_treatment,
    mu,
    n1,
    n1,
    binary = c(5, 6, 7, 8),
    bin_con_resp = bin_con_resp,
    bin_treat_resp = bin_treat_resp
  )

  data_gen_2 <- get_data_gen_2(
    corr_control,
    corr_treatment,
    mu,
    binary = c(5, 6, 7, 8),
    bin_con_resp = bin_con_resp,
    bin_treat_resp = bin_treat_resp
  )

  c(data_gen_1, data_gen_2)
}

get_example_adaption <- function(futility, alt_drop = FALSE) {
  function(design) {
    p <- design$p_values_interim[1:4]

    trt_rej <- design$rej_interim[1:4] & design$rej_interim[5:8]
    if (futility > 0) {
      drop_hyp <- (p > futility) | trt_rej
    } else {
      # drop all but the best performing one (and that as well if it is rejected)
      drop_hyp <- (p != min(p)) | trt_rej
    }

    drop_arms <- which(drop_hyp)

    new_n <- redistribute_n(
      n = c(design$n_treatments[1:4], design$n_controls[1]),
      t = design$t,
      which_drop = drop_arms,
    )

    design |>
      multiarm_drop_arms(
        drop_arms,
        n_cont_2 = rep(new_n[5], 2),
        n_treat_2 = rep(new_n[1:4], 2),
        alt_adj = alt_drop
      )
  }
}


run_example_trial <- function(
  design,
  runs1 = 10,
  runs2 = 100,
  n1 = 50,
  n2 = 50,
  corr = 0.8,
  eff = c(0, 0, 0, 0),
  futility = 0.75,
  alt_drop = TRUE
) {
  dat <- example_data_gen(corr, eff, n1)
  data_gen_1 <- dat[[1]]
  data_gen_2 <- dat[[2]]
  adapt_rule <- get_example_adaption(futility, alt_drop = alt_drop)
  sim_trial(design, runs1, runs2, adapt_rule, data_gen_1, data_gen_2)
}

run_example_trial_bin <- function(
  design,
  runs1 = 10,
  runs2 = 100,
  n1 = 50,
  n2 = 50,
  corr = 0.8,
  eff = c(0, 0, 0, 0),
  futility = 0.75,
  alt_drop = TRUE,
  bin_con_resp,
  bin_treat_resp
) {
  dat <- example_data_gen_bin(corr, eff, n1, bin_con_resp, bin_treat_resp)
  data_gen_1 <- dat[[1]]
  data_gen_2 <- dat[[2]]
  adapt_rule <- get_example_adaption(futility, alt_drop = alt_drop)
  sim_trial(design, runs1, runs2, adapt_rule, data_gen_1, data_gen_2)
}
