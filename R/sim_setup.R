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
    alpha_spending_f = as,
    seq_bonf = FALSE
  )

  design
}

get_sim_data_gen <- function(
  corr,
  eff,
  n1,
  bin,
  bin_con_resp,
  bin_treat_resp
) {
  #the ifelse is not technically necessary, but is more explicit
  if (length(bin) > 0) {
    mu <- c(eff, NA, NA, NA, NA)
  } else {
    mu <- c(eff, eff)
  }

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
    binary = bin,
    bin_con_resp = bin_con_resp,
    bin_treat_resp = bin_treat_resp
  )

  data_gen_2 <- get_data_gen_2(
    corr_control,
    corr_treatment,
    mu,
    binary = bin,
    bin_con_resp = bin_con_resp,
    bin_treat_resp = bin_treat_resp
  )

  c(data_gen_1, data_gen_2)
}

get_sim_adaption <- function(futility, alt_drop = FALSE) {
  function(design) {
    # Only primary p-values are used for decisions
    p <- design$p_values_interim[1:4]

    # arms count as rejected only when primary and secondary hypothesis are
    # rejected
    trt_rej <- design$rej_interim[1:4] & design$rej_interim[5:8]
    if (futility > 0) {
      #only keep arm with primary p above futility
      drop_hyp <- (p > futility) | trt_rej
    } else {
      # drop all but the best performing arm (and that as well if it is rejected)
      drop_hyp <- (p != min(p)) | trt_rej
    }

    drop_arms <- which(drop_hyp)

    #we want the sample size to be redistributed equally across all remaining
    #arms (including control)
    new_n <- redistribute_n(
      n = c(design$n_controls[1], design$n_treatments[1:4]),
      t = design$t,
      which_drop = drop_arms + 1,
    )

    design |>
      multiarm_drop_arms(
        c(drop_arms, drop_arms + 4),
        n_cont_2 = rep(new_n[1], 2),
        n_treat_2 = rep(new_n[2:5], 2),
        alt_adj = alt_drop
      )
  }
}
