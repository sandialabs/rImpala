#' @title Hierarchical Calibration
#'
#' @description Runs hierarchical Bayesian model calibration with adaptive
#'   MCMC, tempering and a joint decorrelation step. Every calibration
#'   parameter vector \eqn{\theta_j} -- one per value of `theta_ind` in each
#'   [addVecExperiments()] call -- gets its own value, and the
#'   \eqn{\theta_j} are tied together by
#'   \eqn{\theta_j \sim N(\theta_0, \Sigma_0)}, truncated to the constraints,
#'   with the hyperpriors set by [setHierPriors()]. Parameters live on the 0-1
#'   scale of `bounds`.
#'
#' @param setup an object of class `CalibSetup`, with [setHierPriors()] called
#'
#' @return a list with the following elements
#' - theta: list over experiments of `nmcmc x ntemps x ntheta x p` arrays of
#'   \eqn{\theta_j} samples (0-1 scale)
#' - theta0: `nmcmc x ntemps x p` array of \eqn{\theta_0} samples
#' - Sigma0: `nmcmc x ntemps x p x p` array of \eqn{\Sigma_0} samples
#' - s2: list over experiments of error variance samples
#' - count: list over experiments of `ntemps x ntheta` acceptance counts
#' - count_s2: acceptance counts for Metropolis error variance updates
#' - count_decor: `ntemps x p` acceptance counts of the decorrelation step
#' - count_temper: `ntemps x ntemps` accepted tempering swaps
#' - cov_theta_cand: adaptive Metropolis state for theta
#' - cov_ls2_cand: adaptive Metropolis state for the log error variances
#' - pred_curr: final predictions for each experiment
#' - llik: log likelihood of the cold chain at each iteration
#' - theta_native: list over experiments of lists over \eqn{j} of data frames
#'   of cold chain \eqn{\theta_j} samples on the native scale
#' - theta0_native: data frame of cold chain \eqn{\theta_0} samples on the
#'   native scale
#' - theta_parent_native: data frame of draws from
#'   \eqn{N(\theta_0, \Sigma_0)} (truncated to the constraints) at each cold
#'   chain sample, i.e. the parameter of a new, unobserved experiment
#'
#' @details Within a temperature the likelihood, the \eqn{N(\theta_0,\Sigma_0)}
#'   density of each \eqn{\theta_j} and the error variance prior are tempered,
#'   while the hyperpriors on \eqn{\theta_0} and \eqn{\Sigma_0} are not; the
#'   Gibbs updates, the decorrelation step and the swaps all target that same
#'   tempered posterior. Discrepancy bases are not supported here.
#'
#'   Error variance groups (`s2_ind`) may be shared across several
#'   \eqn{\theta_j}. Parameters fixed with [fixTheta()] are not supported, and
#'   priors from [addThetaPrior()] are ignored in favour of the hierarchical
#'   prior.
#'
#' @export
#'
calibHier <- function(setup) {
  check_hier_setup(setup, "calibHier")

  nexp = setup$nexp
  ntemps = setup$ntemps
  nmcmc = setup$nmcmc
  p = setup$p
  itl = setup$itl
  ntheta_vec = setup$ntheta
  ntheta_tot = sum(ntheta_vec)

  for (i in 1:nexp) {
    if (setup$models[[i]]$nd > 0) {
      cli::cli_abort("Discrepancy (D) is not supported by {.fn calibHier}; experiment {i} has one.")
    }
  }

  # ---- per-experiment bookkeeping ----
  wt = lapply(1:nexp, function(i) get_wt(setup, i))
  s2_ind_mat = lapply(1:nexp, function(i) outer(setup$s2_ind[[i]], 1:setup$ns2[[i]], "=="))
  ny_w = lapply(1:nexp, function(i) as.numeric(wt[[i]] %*% s2_ind_mat[[i]]))
  which_theta = lapply(1:nexp, function(i) {
    lapply(1:ntheta_vec[i], function(j) which(setup$theta_ind[[i]] == j))
  })

  # ---- hyperpriors ----
  theta0_prior_mean = setup$theta0_prior_mean
  theta0_prior_prec = chol2inv(chol(setup$theta0_prior_cov))
  theta0_prior_ldet = 2 * sum(log(diag(chol(setup$theta0_prior_cov))))
  Sigma0_prior_df = setup$Sigma0_prior_df
  Sigma0_prior_scale = setup$Sigma0_prior_scale

  # ---- storage ----
  theta0 = array(0, c(nmcmc, ntemps, p))
  Sigma0 = array(0, c(nmcmc, ntemps, p, p))
  theta = lapply(ntheta_vec, function(n) array(0, c(nmcmc, ntemps, n, p)))
  log_s2 = lapply(1:nexp, function(i) array(0, c(nmcmc, ntemps, setup$ns2[[i]])))

  theta0[1, , ] = theta0_start_values(setup)
  Sigma0_start = Sigma0_prior_scale / (Sigma0_prior_df - p - 1)
  Sigma0_inv = array(0, c(ntemps, p, p))
  Sigma0_ldet = rep(0, ntemps)
  for (t in 1:ntemps) {
    Sigma0[1, t, , ] = Sigma0_start
    Sigma0_inv[t, , ] = chol2inv(chol(Sigma0_start))
    Sigma0_ldet[t] = 2 * sum(log(diag(chol(Sigma0_start))))
  }

  # current state helpers
  th_now = function(i, m) array(theta[[i]][m, , , ], c(ntemps, ntheta_vec[i], p))
  t0_now = function(m) matrix(theta0[m, , ], ntemps, p)

  # covariance of the observations that belong to theta j, at temperature t
  cov_ij = function(i, ls2_t, j) {
    obs = which_theta[[i]][[j]]
    s2_obs = exp(ls2_t)[setup$s2_ind[[i]][obs]]
    lik_cov_inv(setup$models[[i]], s2_obs, wt = wt[[i]][obs], inds = obs)
  }
  llik_ij = function(i, pred_t, cov, j) {
    obs = which_theta[[i]][[j]]
    as.numeric(llik(setup$models[[i]], setup$ys[[i]][obs], pred_t[obs], cov, wt = wt[[i]][obs]))
  }
  # (ntemps x ntheta) log N(theta_ij; theta0_t, Sigma0_t)
  lprior_theta = function(th, t0) {
    out = matrix(0, ntemps, dim(th)[2])
    for (t in 1:ntemps) {
      out[t, ] = dmvn_rows(matrix(th[t, , ], ncol = p), t0[t, ], Sigma0_inv[t, , ], Sigma0_ldet[t])
    }
    out
  }

  pred_curr = vector("list", nexp)
  llik_curr = vector("list", nexp)
  cov_curr = vector("list", nexp)
  for (i in 1:nexp) {
    log_s2[[i]][1, , ] = matrix(log(setup$sd_est[[i]]^2), ntemps, setup$ns2[[i]], byrow = TRUE)
    init = array(0, c(ntemps, ntheta_vec[i], p))
    for (t in 1:ntemps) {
      L = t(chol(Sigma0[1, t, , ]))
      init[t, , ] = rmvn_constrained(matrix(theta0[1, t, ], ntheta_vec[i], p, byrow = TRUE), L, setup)
    }
    theta[[i]][1, , , ] = init
    pred_curr[[i]] = eval_by_theta(setup$models[[i]], init, setup$theta_ind[[i]], setup)
    llik_curr[[i]] = matrix(0, ntemps, ntheta_vec[i])
    cov_curr[[i]] = vector("list", ntemps)
    for (t in 1:ntemps) {
      cov_curr[[i]][[t]] = vector("list", ntheta_vec[i])
      for (j in 1:ntheta_vec[i]) {
        cov_curr[[i]][[t]][[j]] = cov_ij(i, log_s2[[i]][1, t, ], j)
        llik_curr[[i]][t, j] = llik_ij(i, pred_curr[[i]][t, ], cov_curr[[i]][[t]][[j]], j)
      }
    }
  }

  cov_theta_cand = am_hist_init(ntemps, ntheta_vec, p, setup$start_var_theta,
                                setup$start_adapt_iter, setup$start_tau_theta)
  cov_ls2_cand = lapply(1:nexp, function(i) {
    AMcov_pool(ntemps, setup$ns2[[i]], setup$start_var_ls2,
               setup$start_adapt_iter, setup$start_tau_ls2)
  })

  count = lapply(ntheta_vec, function(n) matrix(0, ntemps, n))
  count_s2 = matrix(0, nexp, ntemps)
  count_decor = matrix(0, ntemps, p)
  count_temper = matrix(0, ntemps, ntemps)
  llik_trace = rep(0, nmcmc)

  pb <- progress::progress_bar$new(
    format = "  MCMC [:bar] :current/:total (:percent) in :eta [:tick_rate it/sec]",
    total = nmcmc,
    clear = FALSE,
    width = 60
  )

  for (m in 2:nmcmc) {
    theta0[m, , ] = theta0[m - 1, , ]
    Sigma0[m, , , ] = Sigma0[m - 1, , , ]

    # ---- emulator step ----
    for (i in 1:nexp) {
      theta[[i]][m, , , ] = theta[[i]][m - 1, , , ]
      log_s2[[i]][m, , ] = log_s2[[i]][m - 1, , ]
      setup$models[[i]] = step_m(setup$models[[i]])
      if (setup$models[[i]]$stochastic) {
        pred_curr[[i]] = eval_by_theta(setup$models[[i]], th_now(i, m), setup$theta_ind[[i]], setup)
        for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
          llik_curr[[i]][t, j] = llik_ij(i, pred_curr[[i]][t, ], cov_curr[[i]][[t]][[j]], j)
        }
      }
    }

    # ---- adaptive Metropolis for every theta_ij ----
    cov_theta_cand = am_hist_update(cov_theta_cand, theta, m)
    cov_theta_cand = am_hist_scale(cov_theta_cand)
    t0_prev = t0_now(m - 1)

    for (i in 1:nexp) {
      cur = th_now(i, m)
      cand = cur
      for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
        cand[t, j, ] = cur[t, j, ] +
          t(chol(cov_theta_cand$S[[i]][t, j, , ])) %*% stats::rnorm(p)
      }
      good = matrix(theta_ok(matrix(cand, ntemps * ntheta_vec[i], p), setup), ntemps, ntheta_vec[i])
      eval_arr = cur
      for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
        if (good[t, j]) eval_arr[t, j, ] = cand[t, j, ]
      }
      pred_cand = eval_by_theta(setup$models[[i]], eval_arr, setup$theta_ind[[i]], setup)
      llik_cand = llik_curr[[i]]
      for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
        if (good[t, j]) {
          llik_cand[t, j] = llik_ij(i, pred_cand[t, ], cov_curr[[i]][[t]][[j]], j)
        }
      }
      alpha = matrix(-Inf, ntemps, ntheta_vec[i])
      lp_diff = lprior_theta(eval_arr, t0_prev) - lprior_theta(cur, t0_prev)
      alpha[good] = (itl * (llik_cand - llik_curr[[i]] + lp_diff))[good]
      accept = log(matrix(stats::runif(length(alpha)), ntemps)) < alpha
      for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
        if (accept[t, j]) {
          theta[[i]][m, t, j, ] = cand[t, j, ]
          obs = which_theta[[i]][[j]]
          pred_curr[[i]][t, obs] = pred_cand[t, obs]
          llik_curr[[i]][t, j] = llik_cand[t, j]
        }
      }
      count[[i]] = count[[i]] + accept
      cov_theta_cand$count_100[[i]] = cov_theta_cand$count_100[[i]] + accept
    }
    cov_theta_cand = am_hist_tau(cov_theta_cand, m)

    # ---- error variances ----
    for (i in 1:nexp) {
      s2mode = setup$models[[i]]$s2
      if (s2mode == "gibbs") {
        dev_sq = weighted_dev_sq(pred_curr[[i]], setup$ys[[i]], wt[[i]], s2_ind_mat[[i]])
        for (t in 1:ntemps) {
          log_s2[[i]][m, t, ] = gibbs_log_s2(dev_sq[t, ], ny_w[[i]], itl[t],
                                             setup$ig_a[[i]], setup$ig_b[[i]],
                                             setup$sd_lower[[i]], setup$sd_upper[[i]])
          for (j in 1:ntheta_vec[i]) {
            cov_curr[[i]][[t]][[j]] = cov_ij(i, log_s2[[i]][m, t, ], j)
            llik_curr[[i]][t, j] = llik_ij(i, pred_curr[[i]][t, ], cov_curr[[i]][[t]][[j]], j)
          }
        }
      } else if (s2mode == "fix") {
        next
      } else {
        cov_ls2_cand[[i]] = update_m(cov_ls2_cand[[i]], log_s2[[i]], m)
        ls2_cand = gen_cand(cov_ls2_cand[[i]], log_s2[[i]], m)
        ls2_curr = matrix(log_s2[[i]][m - 1, , ], ntemps, setup$ns2[[i]])
        llik_candi = matrix(0, ntemps, ntheta_vec[i])
        cov_candi = vector("list", ntemps)
        for (t in 1:ntemps) {
          cov_candi[[t]] = vector("list", ntheta_vec[i])
          ls2_t = pmin(ls2_cand[t, ], log(1e100))
          for (j in 1:ntheta_vec[i]) {
            cov_candi[[t]][[j]] = cov_ij(i, ls2_t, j)
            llik_candi[t, j] = llik_ij(i, pred_curr[[i]][t, ], cov_candi[[t]][[j]], j)
          }
        }
        alpha_s2 = itl * (rowSums(llik_candi) - rowSums(llik_curr[[i]]))
        alpha_s2 = alpha_s2 + itl * s2_kern_sum(setup$s2_prior_kern[[i]], ls2_cand, setup$ig_a[[i]], setup$ig_b[[i]])
        alpha_s2 = alpha_s2 + itl * ls2_rowsum(ls2_cand)
        alpha_s2 = alpha_s2 - itl * s2_kern_sum(setup$s2_prior_kern[[i]], ls2_curr, setup$ig_a[[i]], setup$ig_b[[i]])
        alpha_s2 = alpha_s2 - itl * ls2_rowsum(ls2_curr)
        alpha_s2[!ls2_in_bounds(ls2_cand, setup$sd_lower[[i]], setup$sd_upper[[i]])] = -Inf
        for (t in which(log(stats::runif(ntemps)) < alpha_s2)) {
          count_s2[i, t] = count_s2[i, t] + 1
          log_s2[[i]][m, t, ] = ls2_cand[t, ]
          llik_curr[[i]][t, ] = llik_candi[t, ]
          cov_curr[[i]][[t]] = cov_candi[[t]]
          cov_ls2_cand[[i]]$count_100[t] = cov_ls2_cand[[i]]$count_100[t] + 1
        }
        cov_ls2_cand[[i]] = update_tau(cov_ls2_cand[[i]], m)
      }
    }

    # ---- Gibbs update theta0 (truncated to the constraints) ----
    theta_sum = matrix(0, ntemps, p)
    for (i in 1:nexp) {
      cur = th_now(i, m)
      for (t in 1:ntemps) theta_sum[t, ] = theta_sum[t, ] + colSums(matrix(cur[t, , ], ncol = p))
    }
    for (t in 1:ntemps) {
      prec = ntheta_tot * itl[t] * Sigma0_inv[t, , ] + theta0_prior_prec
      cc = chol2inv(chol(prec))
      dd = itl[t] * Sigma0_inv[t, , ] %*% theta_sum[t, ] + theta0_prior_prec %*% theta0_prior_mean
      theta0[m, t, ] = rmvn_constrained(t(cc %*% dd), t(chol(cc)), setup)
    }

    # ---- Gibbs update Sigma0 ----
    for (t in 1:ntemps) {
      M = matrix(0, p, p)
      for (i in 1:nexp) {
        d = sweep(matrix(theta[[i]][m, t, , ], ncol = p), 2, theta0[m, t, ])
        M = M + crossprod(d)
      }
      S_t = rinvwishart(Sigma0_prior_df + ntheta_tot * itl[t], Sigma0_prior_scale + itl[t] * M)
      Sigma0[m, t, , ] = S_t
      R = chol(S_t)
      Sigma0_inv[t, , ] = chol2inv(R)
      Sigma0_ldet[t] = 2 * sum(log(diag(R)))
    }

    # ---- joint decorrelation: shift theta0 and every theta_j together ----
    if (m %% setup$decor == 0) {
      for (k in 1:p) {
        z = stats::rnorm(1) * 0.1
        t0_cur = t0_now(m)
        t0_cand = t0_cur
        t0_cand[, k] = t0_cand[, k] + z
        alpha_tot = ifelse(theta_ok(t0_cand, setup), 0, -Inf)
        alpha_tot = alpha_tot -
          0.5 * rowSums((sweep(t0_cand, 2, theta0_prior_mean) %*% theta0_prior_prec) * sweep(t0_cand, 2, theta0_prior_mean)) +
          0.5 * rowSums((sweep(t0_cur, 2, theta0_prior_mean) %*% theta0_prior_prec) * sweep(t0_cur, 2, theta0_prior_mean))
        cands = vector("list", nexp)
        pred_cands = vector("list", nexp)
        llik_cands = vector("list", nexp)
        for (i in 1:nexp) {
          cur = th_now(i, m)
          cand = cur
          cand[, , k] = cand[, , k] + z
          good = matrix(theta_ok(matrix(cand, ntemps * ntheta_vec[i], p), setup), ntemps, ntheta_vec[i])
          alpha_tot[rowSums(!good) > 0] = -Inf
          eval_arr = cur
          for (t in 1:ntemps) if (all(good[t, ])) eval_arr[t, , ] = cand[t, , ]
          pred_cands[[i]] = eval_by_theta(setup$models[[i]], eval_arr, setup$theta_ind[[i]], setup)
          llik_cands[[i]] = llik_curr[[i]]
          for (t in 1:ntemps) {
            if (all(good[t, ])) {
              for (j in 1:ntheta_vec[i]) {
                llik_cands[[i]][t, j] = llik_ij(i, pred_cands[[i]][t, ], cov_curr[[i]][[t]][[j]], j)
              }
            }
          }
          lp_diff = lprior_theta(cand, t0_cand) - lprior_theta(cur, t0_cur)
          alpha_tot = alpha_tot + itl * rowSums(llik_cands[[i]] - llik_curr[[i]] + lp_diff)
          cands[[i]] = cand
        }
        accept = which(log(stats::runif(ntemps)) < alpha_tot)
        for (t in accept) {
          theta0[m, t, ] = t0_cand[t, ]
          for (i in 1:nexp) {
            theta[[i]][m, t, , ] = cands[[i]][t, , ]
            pred_curr[[i]][t, ] = pred_cands[[i]][t, ]
            llik_curr[[i]][t, ] = llik_cands[[i]][t, ]
          }
          count_decor[t, k] = count_decor[t, k] + 1
        }
      }
    }

    # ---- tempering swaps ----
    if ((m > setup$start_temper) && (ntemps > 1)) {
      t0_cur = t0_now(m)
      for (sweep_k in 1:setup$nswap) {
        sw = matrix(sample(1:ntemps, 2 * setup$nswap_per, FALSE), 2)
        # tempered part of the log posterior at each temperature's current state
        tempered = rep(0, ntemps)
        for (i in 1:nexp) {
          ls2 = matrix(log_s2[[i]][m, , ], ntemps, setup$ns2[[i]])
          tempered = tempered + rowSums(llik_curr[[i]]) +
            s2_kern_sum(setup$s2_prior_kern[[i]], ls2, setup$ig_a[[i]], setup$ig_b[[i]]) +
            rowSums(lprior_theta(th_now(i, m), t0_cur))
        }
        sw_alpha = (itl[sw[2, ]] - itl[sw[1, ]]) * (tempered[sw[1, ]] - tempered[sw[2, ]])
        for (s in which(log(stats::runif(setup$nswap_per)) < sw_alpha)) {
          a = sw[1, s]
          b = sw[2, s]
          count_temper[a, b] = count_temper[a, b] + 1
          for (i in 1:nexp) {
            theta[[i]][m, c(a, b), , ] = theta[[i]][m, c(b, a), , ]
            log_s2[[i]][m, c(a, b), ] = log_s2[[i]][m, c(b, a), ]
            pred_curr[[i]] = swap_rows(pred_curr[[i]], a, b)
            llik_curr[[i]] = swap_rows(llik_curr[[i]], a, b)
            cov_curr[[i]][c(a, b)] = cov_curr[[i]][c(b, a)]
          }
          theta0[m, c(a, b), ] = theta0[m, c(b, a), ]
          Sigma0[m, c(a, b), , ] = Sigma0[m, c(b, a), , ]
          Sigma0_inv = swap_rows(Sigma0_inv, a, b)
          Sigma0_ldet = swap_rows(Sigma0_ldet, a, b)
          t0_cur = t0_now(m)
        }
      }
    }

    llik_trace[m] = sum(vapply(llik_curr, function(x) sum(x[1, ]), numeric(1)))
    pb$tick()
  }

  # ---- outputs ----
  s2 = lapply(log_s2, exp)
  pn = names(setup$bounds)
  to_native = function(x) {
    out = data.frame(tran_unif(matrix(x, ncol = p), setup$bounds_mat, pn))
    names(out) = pn
    out
  }
  theta_native = lapply(1:nexp, function(i) {
    lapply(1:ntheta_vec[i], function(j) to_native(theta[[i]][, 1, j, ]))
  })
  theta0_native = to_native(theta0[, 1, ])
  parent = matrix(0, nmcmc, p)
  for (mm in 1:nmcmc) {
    parent[mm, ] = rmvn_constrained(theta0[mm, 1, ], t(chol(Sigma0[mm, 1, , ])), setup)
  }

  count_temper = count_temper + t(count_temper) - diag(diag(count_temper), ntemps)

  list(
    theta = theta,
    theta0 = theta0,
    Sigma0 = Sigma0,
    s2 = s2,
    count = count,
    count_s2 = count_s2,
    count_decor = count_decor,
    count_temper = count_temper,
    cov_theta_cand = cov_theta_cand,
    cov_ls2_cand = cov_ls2_cand,
    pred_curr = pred_curr,
    llik = llik_trace,
    theta_native = theta_native,
    theta0_native = theta0_native,
    theta_parent_native = to_native(parent)
  )
}
