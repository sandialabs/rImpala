#' This function runs pooled Bayesian Model Calibration with adaptive MCMC,
#' tempering, and decorrelation steps
#'
#' input theta will be normalized to 0-1 and sampled from uniform priors
#'
#' @param setup: an object of class `CalibSetup`
#'
#' @return a list with the following elements
#'
#' - theta: mcmc samples of variables
#' - s2: mcmc samples of error variance
#' - count: number of counts of acceptance
#' - count_s2: number of counts of acceptance on error variance
#' - count_decor: number of times decorrelation occurred
#' - cov_theta_cand: final theta covariance
#' - cov_ls2_cand: final error covariance
#' - pred_curr: current emulator predictions
#' - discrep_vars: discrepancy coefficients
#' - llik: log likelihood
#' - theta_native: mcmc samples of variables in native scale
#'
#' @export
#'
calibPool <- function(setup) {
  theta = array(0, dim = c(setup$nmcmc, setup$ntemps, setup$p))
  log_s2 <- vector(mode = "list", length = setup$nexp)
  s2_ind_mat <- vector(mode = "list", length = setup$nexp)

  for (i in 1:setup$nexp) {
    log_s2[[i]] = array(1, dim = c(setup$nmcmc, setup$ntemps, setup$ns2[[i]]))
    s2_ind_mat[[i]] = setup$s2_ind[[i]] == 1:setup$ns2[[i]]
  }

  theta_start = matrix(runif(setup$ntemps * setup$p), setup$ntemps)

  good = setup$checkConstraints(tran_unif(theta_start, setup$bounds_mat, names(setup$bounds)),
                                setup$bounds)

  while (any(!good)) {
    theta_start[!good, ] = matrix(runif(sum(!good) * setup$p), sum(!good))
    good[!good] = setup$checkConstraints(tran_unif(theta_start[!good, ], setup$bounds_mat, names(setup$bounds)),
                                         setup$bounds)
  }

  theta[1, , ] = theta_start

  # matrix of temperatures for use with alpha calculation--to skip nested for loops.
  itl_mat = vector(mode = "list", length = setup$nexp)
  for (i in 1:setup$nexp) {
    itl_mat[[i]] = (matrix(1, setup$ns2[[i]], setup$ntemps) * setup$itl)
  }

  pred_curr = vector(mode = "list", length = setup$nexp)
  llik_curr = matrix(0, setup$nexp, setup$ntemps)
  marg_lik_cov_cur = vector(mode = "list", length = setup$nexp)
  for (i in 1:setup$nexp) {
    marg_lik_cov_cur[[i]] = vector(mode = "list", length = setup$ntemps)
    for (t in 1:setup$ntemps) {
      tmp = exp(log_s2[[i]][1, t, setup$s2_ind[[i]]])
      marg_lik_cov_cur[[i]][[t]] = lik_cov_inv(setup$models[[i]], tmp[setup$s2_ind[[i]]])
    }
  }

  for (i in 1:setup$nexp) {
    pred_curr[[i]] = evalm(setup$models[[i]],
                          tran_unif(theta[1, , ], setup$bounds_mat, names(setup$bounds)),
                          TRUE)

    for (t in 1:setup$ntemps) {
      llik_curr[i, t] = llik(setup$models[[i]],
                             setup$ys[[i]],
                             pred_curr[[i]][t, ],
                             marg_lik_cov_cur[[i]][[t]])
    }
  }

  cov_theta_cand = AMcov_pool(
    setup$ntemps,
    setup$p,
    setup$start_var_theta,
    setup$start_adapt_iter,
    setup$start_tau_theta
  )
  cov_ls2_cand = vector(mode = "list", length = setup$nexp)
  for (i in 1:setup$nexp) {
    cov_ls2_cand[[i]] = AMcov_pool(
      setup$ntemps,
      setup$ns2[[i]],
      setup$start_var_ls2,
      setup$start_adapt_iter,
      setup$start_tau_ls2
    )
  }

  count = matrix(0, setup$ntemps, setup$ntemps)
  count_s2 = matrix(0, setup$nexp, setup$ntemps)
  count_decor = matrix(0, setup$p, setup$ntemps)

  discrep_curr = pred_curr
  for (i in 1:length(discrep_curr)) {
    discrep_curr[[i]] = discrep_curr[[i]] * 0
  }

  discrep_vars = vector(mode = "list", length = setup$nexp)
  for (i in 1:setup$nexp) {
    discrep_vars[[i]] = array(0, dim = c(setup$nmcmc, setup$ntemps, setup$models[[i]]$nd))
  }

  alpha = rep(1, setup$ntemps) * -Inf
  sw_alpha = rep(0, setup$nswap_per)

  llik = rep(0, setup$nmcmc)

  # start MCMC
  pb <- progress::progress_bar$new(
    format = "  MCMC [:bar] :current/:total (:percent) in :eta [:tick_rate it/sec]",
    total = setup$nmcmc,
    clear = FALSE,
    width = 60
  )
  for (m in 2:setup$nmcmc) {
    theta[m, , ] = theta[m - 1, , ]

    for (i in 1:setup$nexp) {
      log_s2[[i]][m, , ] = log_s2[[i]][m - 1, , ]
      if (setup$models[[i]]$nd > 0) {
        for (t in 1:setup$ntemps) {
          discrep_vars[[i]][m, t, ] = discrep_sample(setup$models[[i]],
                                                     setup$ys[[i]],
                                                     pred_curr[[i]][t, ],
                                                     marg_lik_cov_cur[[i]][[t]],
                                                     setup$itl[t])
          discrep_curr[[i]][t, ] = setup$models[[i]]$D %*% discrep_vars[[i]][m, t, ]
        }
      }

      setup$models[[i]] = step_m(setup$models[[i]])
      if (setup$models[[i]]$stochastic) {
        pred_curr[[i]] = evalm(setup$models[[i]],
                              tran_unif(theta[m, , ], setup$bounds_mat, names(setup$bounds)),
                              TRUE)
      }
      if ((setup$models[[i]]$nd > 0) |
          (setup$models[[i]]$stochastic)) {
        for (t in 1:setup$ntemps) {
          llik_curr[i, t] = llik(
            setup$models[[i]],
            setup$ys[[i]] - discrep_curr[[i]][t, ],
            pred_curr[[i]][t, ],
            marg_lik_cov_cur[[i]][[t]]
          )
        }
      }
    }

    # adaptive Metropolis for each temperature
    cov_theta_cand = update_m(cov_theta_cand, theta, m)

    # generate proposal
    theta_cand = gen_cand(cov_theta_cand, theta, m)
    good_values = setup$checkConstraints(tran_unif(theta_cand, setup$bounds_mat, names(setup$bounds)),
                                         setup$bounds)

    # get predictions and SSE
    pred_cand = pred_curr
    llik_cand = llik_curr
    if (any(good_values)) {
      llik_cand[, good_values = 0]
      for (i in 1:setup$nexp) {
        theta_tmp = matrix(theta_cand[good_values, ], ncol = setup$p)
        pred_cand[[i]][good_values, ] = evalm(setup$models[[i]],
                                             tran_unif(theta_tmp, setup$bounds_mat, names(setup$bounds)),
                                             TRUE)
        for (t in 1:setup$ntemps) {
          llik_cand[i, t] = llik(
            setup$models[[i]],
            setup$ys[[i]] - discrep_curr[[i]][t, ],
            pred_cand[[i]][t, ],
            marg_lik_cov_cur[[i]][[t]]
          )
        }
      }
    }

    llik_diff = (colSums(llik_cand) - colSums(llik_curr))
    llik_diff = llik_diff[good_values]

    alpha = rep(1, setup$ntemps) * -Inf
    alpha[good_values] = setup$itl[good_values] * llik_diff
    idx = which(log(runif(setup$ntemps)) < alpha)
    for (t in idx) {
      theta[m, t, ] = theta_cand[t, ]
      count[t, t] = count[t, t] + 1
      for (i in 1:setup$nexp) {
        llik_curr[i, t] = llik_cand[i, t]
        pred_curr[[i]][t, ] = pred_cand[[i]][t, ]
      }
      cov_theta_cand$count_100[t] = cov_theta_cand$count_100[t] + 1
    }

    # diminishing adaptation based on acceptance rate for each temperature
    cov_theta_cand = update_tau(cov_theta_cand, m)

    # decorrelation step
    if (m %% setup$decor == 0) {
      for (k in 1:setup$p) {
        theta_cand = theta[m, , ]
        theta_cand[, k] = runif(setup$ntemps)
        good_values = setup$checkConstraints(tran_unif(theta_cand, setup$bounds_mat, names(setup$bounds)),
                                             setup$bounds)
        pred_cand = pred_curr
        llik_cand = llik_curr

        if (any(good_values)) {
          llik_cand[, good_values] = 0
          for (i in 1:setup$nexp) {
            theta_tmp = matrix(theta_cand[good_values, ], ncol = setup$p)
            pred_cand[[i]][good_values, ] = evalm(setup$models[[i]],
                                                 tran_unif(theta_tmp, setup$bounds_mat, names(setup$bounds)),
                                                 TRUE)
            for (t in 1:setup$ntemps) {
              llik_cand[i, t] = llik(
                setup$models[[i]],
                setup$ys[[i]] - discrep_curr[[i]][t, ],
                pred_cand[[i]][t, ],
                marg_lik_cov_cur[[i]][[t]]
              )

            }
          }
        }

        alpha = rep(1, setup$ntemps) * -Inf

        llik_diff = (colSums(llik_cand) - colSums(llik_curr))
        llik_diff = llik_diff[good_values]

        alpha[good_values] = setup$itl[good_values] * llik_diff

        idx = which(log(runif(setup$ntemps)) < alpha)
        for (t in idx) {
          theta[m, t, k] = theta_cand[t, k]
          count_decor[k, t] = count_decor[k, t] + 1
          for (i in 1:setup$nexp) {
            pred_curr[[i]][t, ] = pred_cand[[i]][t, ]
            llik_curr[i, t] = llik_cand[i, t]
          }
        }
      }
    }

    # update s2
    for (i in 1:setup$nexp) {
      if (setup$models[[i]]$s2 == 'gibbs') {
        # gibbs update s2
        dev_sq = (pred_curr[[i]] - setup$ys[[i]]) ^ 2 %*% s2_ind_mat[[i]]
        log_s2[[i]][m, ] = log(1 / rgamma(
          itl_mat[[i]] * (setup$ny_s2[[i]] / 2 + setup$ig_a[[i]] + 1) - 1,
          1 / (itl_mat[[i]] * (setup$ig_b[[i]] + dev_sq / 2))
        ))
        for (t in 1:setup$ntemps) {
          tmpi = exp(log_s2[[i]][m, t])
          marg_lik_cov_cur[[i]][[t]] = lik_cov_inv(setup$models[[i]], tmpi[setup$s2_ind[i]])
          llik_curr[i, t] = llik(
            setup$models[[i]],
            setup$ys[[i]] - discrep_curr[[i]][t, ],
            pred_curr[[i]][t, ],
            marg_lik_cov_cur[[i]][[t]]
          )
        }
      } else {
        # M-H update s2
        cov_ls2_cand[[i]] = update_m(cov_ls2_cand[[i]], log_s2[[i]], m)
        ls2_candi = gen_cand(cov_ls2_cand[[i]], log_s2[[i]], m)

        llik_candi = rep(0, setup$ntemps)
        marg_lik_cov_candi = vector(mode = "list", length = setup$ntemps)

        for (t in 1:setup$ntemps) {
          tmpi = exp(ls2_candi[t])
          marg_lik_cov_candi[[t]] = lik_cov_inv(setup$models[[i]], tmpi[setup$s2_ind[[i]]])
          llik_candi[t] = llik(
            setup$models[[i]],
            setup$ys[[i]] - discrep_curr[[i]][t, ],
            pred_curr[[i]][t, ],
            marg_lik_cov_candi[[t]]
          )
        }

        llik_diffi = (llik_candi - llik_curr[i, ])
        alpha_s2 = setup$itl * llik_diffi
        alpha_s2 = alpha_s2 + setup$itl * rowSums(setup$s2_prior_kern[[i]](exp(ls2_candi), setup$ig_a[[i]], setup$ig_b[[i]]))
        alpha_s2 = alpha_s2 + setup$itl * rowSums(ls2_candi)
        alpha_s2 = alpha_s2 - setup$itl * colSums(setup$s2_prior_kern[[i]](exp(t(log_s2[[i]][m - 1, , ])), setup$ig_a[[i]], setup$ig_b[[i]]))
        alpha_s2 = alpha_s2 - setup$itl * colSums(t(log_s2[[i]][m - 1, , ]))

        idx = which(log(runif(setup$ntemps)) < alpha_s2)
        for (t in idx) {
          count_s2[i, t] = count_s2[i, t] + 1
          llik_curr[i, t] = llik_candi[t]
          log_s2[[i]][m, t, ] = ls2_candi[t]
          marg_lik_cov_cur[[i]][[t]] = marg_lik_cov_candi[[t]]
          cov_ls2_cand[[i]]$count_100 = cov_ls2_cand[[i]]$count_100 + 1
        }

        cov_ls2_cand[[i]] = update_tau(cov_ls2_cand[[i]], m)
      }
    }

    # tempering swaps
    if ((m > setup$start_temper) & (setup$ntemps > 1)) {
      for (k in 1:setup$nswap) {
        sw = sample(1:setup$ntemps, 2 * setup$nswap_per, FALSE)
        sw = matrix(sw, setup$nswap_per, 2)
        sw = t(sw)

        sw_alpha = rep(0, setup$nswap_per)
        sw_alpha = sw_alpha + (setup$itl[sw[2, ]] - setup$itl[sw[1, ]]) * (colSums(llik_curr[, sw[1, ]]) - colSums(llik_curr[, sw[2, ]]))
        for (i in 1:setup$nexp) {
          sw_alpha = sw_alpha + (setup$itl[sw[2, ]] - setup$itl[sw[1, ]]) *
            (rowSums(matrix(
              setup$s2_prior_kern[[i]](exp(log_s2[[i]][m, sw[1, ], ]), setup$ig_a[[i]], setup$ig_b[[i]])
            )) -
              rowSums(matrix(
                setup$s2_prior_kern[[i]](exp(log_s2[[i]][m, sw[2, ], ]), setup$ig_a[[i]], setup$ig_b[[i]])
              )))
          if (setup$models[[i]]$nd > 0) {
            sw_alpha = sw_alpha + (setup$itl[sw[2, ]] - setup$itl[sw[1, ]]) *
              (
                -0.5 * rowSums(discrep_vars[[i]][m, (sw[1, ])] ^ 2) / setup$models[[i]]$discrep_tau +
                  0.5 * rowSums(discrep_vars[[i]][m, (sw[2, ])] ^
                                  2) / setup$models[[i]]$discrep_tau
              )
          }
        }

        idx = which(log(runif(setup$nswap_per)) < sw_alpha)

        for (tti in idx) {
          tt = sw[, tti]
          for (i in 1:setup$nexp) {
            log_s2[[i]][m, tt[1], ] = log_s2[[i]][m, tt[2], ]
            log_s2[[i]][m, tt[2], ] = log_s2[[i]][m, tt[1], ]
            marg_lik_cov_cur[[i]][tt[1]] = marg_lik_cov_cur[[i]][tt[2]]
            marg_lik_cov_cur[[i]][tt[2]] = marg_lik_cov_cur[[i]][tt[1]]
            pred_curr[[i]][tt[1], ] = pred_curr[[i]][tt[2], ]
            pred_curr[[i]][tt[2], ] = pred_curr[[i]][tt[1], ]
            if (setup$models[[i]]$nd > 0) {
              discrep_curr[[i]][tt[1, ]] = discrep_curr[[i]][tt[2], ]
              discrep_curr[[i]][tt[2, ]] = discrep_curr[[i]][tt[1], ]
              discrep_vars[[i]][m, tt[1]] = discrep_vars[[i]][m, tt[2]]
              discrep_vars[[i]][m, tt[2]] = discrep_vars[[i]][m, tt[1]]
            }
            llik_curr[i, tt[1]] = llik_curr[i, tt[2]]
            llik_curr[i, tt[2]] = llik_curr[i, tt[1]]
          }
          count[tt[1], tt[2]] = count[tt[1], tt[2]] + 1
          theta[m, tt[1], ] = theta[m, tt[2], ]
          theta[m, tt[2], ] = theta[m, tt[1], ]
        }
      }
    }

    llik[m] = sum(llik_curr[, 1])

    pb$tick()
  }

  s2 = log_s2
  for (i in 1:setup$nexp) {
    s2[[i]] = exp(log_s2[[i]])
  }

  theta_native = tran_unif(theta[, 1, ], setup$bounds_mat, names(setup$bounds))

  out <- list(
    theta = theta,
    s2 = s2,
    count = count,
    count_s2 = count_s2,
    count_decor = count_decor,
    cov_theta_cand = cov_theta_cand,
    cov_ls2_cand = cov_ls2_cand,
    pred_curr = pred_curr,
    discrep_vars = discrep_vars,
    llik = llik,
    theta_native = data.frame(theta_native)
  )

  out
}
