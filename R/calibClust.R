#' @title Clustered Calibration
#'
#' @description Runs clustered Bayesian model calibration: the calibration
#'   parameter vectors \eqn{\theta_j} -- one per value of `theta_ind` in each
#'   [addVecExperiments()] call -- are grouped by a Dirichlet process into at
#'   most `nclustmax` clusters that share one parameter value. Cluster values
#'   are drawn from \eqn{N(\theta_0, \Sigma_0)} truncated to the constraints,
#'   with the hyperpriors of [setHierPriors()], and the concentration
#'   \eqn{\eta} has the Gamma prior of [setClusterPriors()]. Parameters live on
#'   the 0-1 scale of `bounds`.
#'
#' @param setup an object of class `CalibSetup`, with [setHierPriors()] and
#'   [setClusterPriors()] called
#'
#' @return a list with the following elements
#' - theta: `nmcmc x ntemps x nclustmax x p` array of cluster values (0-1
#'   scale); only clusters with members are meaningful
#' - theta_hist: list over experiments of `nmcmc x ntemps x ntheta x p` arrays
#'   holding the value each \eqn{\theta_j} took through its cluster
#' - delta: list over experiments of `nmcmc x ntemps x ntheta` arrays of
#'   cluster memberships (1-based)
#' - eta: `nmcmc x ntemps` matrix of the concentration parameter
#' - theta0, Sigma0: hyperparameter samples, as in [calibHier()]
#' - s2: list over experiments of error variance samples
#' - count: accepted cluster value updates per temperature
#' - count_s2: acceptance counts for Metropolis error variance updates
#' - count_temper: `ntemps x ntemps` accepted tempering swaps
#' - pred_curr: final predictions for each experiment
#' - llik: log likelihood of the cold chain at each iteration
#' - nclustmax: the maximum number of clusters
#' - theta_native: list over experiments of lists over \eqn{j} of data frames
#'   of cold chain \eqn{\theta_j} samples on the native scale
#' - theta0_native, theta_parent_native: as in [calibHier()]
#'
#' @details Memberships are updated one \eqn{\theta_j} at a time with the
#'   Chinese restaurant process weights, using the empty clusters (whose values
#'   are redrawn from \eqn{N(\theta_0,\Sigma_0)} every iteration) as the
#'   candidate new clusters; \eqn{\eta} is updated with the auxiliary variable
#'   scheme of Escobar and West (1995). Within a temperature the likelihood,
#'   the error variance prior and the \eqn{N(\theta_0,\Sigma_0)} density of the
#'   occupied clusters are tempered; the cluster process and the hyperpriors
#'   are not.
#'
#'   Discrepancy bases and [fixTheta()] are not supported; priors from
#'   [addThetaPrior()] are ignored.
#'
#' @export
#'
calibClust <- function(setup) {
  check_hier_setup(setup, "calibClust")

  nexp = setup$nexp
  ntemps = setup$ntemps
  nmcmc = setup$nmcmc
  p = setup$p
  itl = setup$itl
  ntheta_vec = setup$ntheta
  ndat = sum(ntheta_vec)
  K = setup$nclustmax
  if (is.null(K)) K = max(ndat, 10)

  for (i in 1:nexp) {
    if (setup$models[[i]]$nd > 0) {
      cli::cli_abort("Discrepancy (D) is not supported by {.fn calibClust}; experiment {i} has one.")
    }
  }

  wt = lapply(1:nexp, function(i) get_wt(setup, i))
  s2_ind_mat = lapply(1:nexp, function(i) outer(setup$s2_ind[[i]], 1:setup$ns2[[i]], "=="))
  ny_w = lapply(1:nexp, function(i) as.numeric(wt[[i]] %*% s2_ind_mat[[i]]))
  which_theta = lapply(1:nexp, function(i) {
    lapply(1:ntheta_vec[i], function(j) which(setup$theta_ind[[i]] == j))
  })
  ny = vapply(setup$ys, length, numeric(1))

  theta0_prior_mean = setup$theta0_prior_mean
  theta0_prior_prec = chol2inv(chol(setup$theta0_prior_cov))
  Sigma0_prior_df = setup$Sigma0_prior_df
  Sigma0_prior_scale = setup$Sigma0_prior_scale
  eta_shape = setup$eta_prior_shape
  eta_rate = setup$eta_prior_rate

  # ---- storage ----
  theta = array(0, c(nmcmc, ntemps, K, p))
  theta0 = array(0, c(nmcmc, ntemps, p))
  Sigma0 = array(0, c(nmcmc, ntemps, p, p))
  eta = matrix(0, nmcmc, ntemps)
  delta = lapply(ntheta_vec, function(n) array(0L, c(nmcmc, ntemps, n)))
  theta_hist = lapply(ntheta_vec, function(n) array(0, c(nmcmc, ntemps, n, p)))
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
  eta[1, ] = eta_shape / eta_rate
  for (t in 1:ntemps) {
    theta[1, t, , ] = rmvn_constrained(matrix(theta0[1, t, ], K, p, byrow = TRUE),
                                       t(chol(Sigma0_start)), setup)
  }
  for (i in 1:nexp) {
    delta[[i]][1, , ] = sample.int(K, ntemps * ntheta_vec[i], replace = TRUE)
    log_s2[[i]][1, , ] = matrix(log(setup$sd_est[[i]]^2), ntemps, setup$ns2[[i]], byrow = TRUE)
  }

  # ---- state helpers ----
  th_now = function(m) array(theta[m, , , ], c(ntemps, K, p))
  dl_now = function(i, m) matrix(delta[[i]][m, , ], ntemps, ntheta_vec[i])
  # (ntemps x ntheta x p) value each theta_j takes through its cluster
  hist_from = function(th, dl) {
    out = array(0, c(ntemps, ncol(dl), p))
    for (t in 1:ntemps) out[t, , ] = th[t, dl[t, ], ]
    out
  }
  extant = function(m) {
    ext = matrix(FALSE, ntemps, K)
    for (i in 1:nexp) {
      dl = dl_now(i, m)
      for (t in 1:ntemps) ext[t, dl[t, ]] = TRUE
    }
    ext
  }
  member_count = function(m) {
    cnt = matrix(0L, ntemps, K)
    for (i in 1:nexp) {
      dl = dl_now(i, m)
      for (t in 1:ntemps) cnt[t, ] = cnt[t, ] + tabulate(dl[t, ], K)
    }
    cnt
  }
  # (ntemps x K x ny) predictions at every cluster value
  eval_clusters = function(i, th) {
    full = evalm(setup$models[[i]],
                 tran_unif(matrix(th, ntemps * K, p), setup$bounds_mat, names(setup$bounds)),
                 TRUE)
    array(matrix(full, ntemps * K, ny[i]), c(ntemps, K, ny[i]))
  }
  cov_ij = function(i, ls2_t, j) {
    obs = which_theta[[i]][[j]]
    s2_obs = exp(ls2_t)[setup$s2_ind[[i]][obs]]
    lik_cov_inv(setup$models[[i]], s2_obs, wt = wt[[i]][obs], inds = obs)
  }
  llik_ij = function(i, pred_t, cov, j) {
    obs = which_theta[[i]][[j]]
    as.numeric(llik(setup$models[[i]], setup$ys[[i]][obs], pred_t[obs], cov, wt = wt[[i]][obs]))
  }
  # prediction of each experiment given memberships: (ntemps x ny)
  pred_from_delta = function(i, pdelta, dl) {
    out = matrix(0, ntemps, ny[i])
    for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
      obs = which_theta[[i]][[j]]
      out[t, obs] = pdelta[t, dl[t, j], obs]
    }
    out
  }
  # llik of every theta_j at every cluster value: (ntemps x K x ntheta)
  llik_all = function(i, pdelta) {
    out = array(0, c(ntemps, K, ntheta_vec[i]))
    for (t in 1:ntemps) for (k in 1:K) for (j in 1:ntheta_vec[i]) {
      out[t, k, j] = llik_ij(i, pdelta[t, k, ], cov_curr[[i]][[t]][[j]], j)
    }
    out
  }
  llik_given_delta = function(i, ldelta, dl) {
    out = matrix(0, ntemps, ntheta_vec[i])
    for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) out[t, j] = ldelta[t, dl[t, j], j]
    out
  }
  lprior_clusters = function(th) {
    out = matrix(0, ntemps, K)
    for (t in 1:ntemps) {
      out[t, ] = dmvn_rows(matrix(th[t, , ], ncol = p), theta0[m_cur, t, ], Sigma0_inv[t, , ], Sigma0_ldet[t])
    }
    out
  }

  # ---- initial predictions and likelihoods ----
  m_cur = 1
  cov_curr = vector("list", nexp)
  pred_delta = vector("list", nexp)
  llik_delta = vector("list", nexp)
  for (i in 1:nexp) {
    theta_hist[[i]][1, , , ] = hist_from(th_now(1), dl_now(i, 1))
    cov_curr[[i]] = lapply(1:ntemps, function(t) {
      lapply(1:ntheta_vec[i], function(j) cov_ij(i, log_s2[[i]][1, t, ], j))
    })
    pred_delta[[i]] = eval_clusters(i, th_now(1))
    llik_delta[[i]] = llik_all(i, pred_delta[[i]])
  }

  am = am_hist_init(ntemps, ntheta_vec, p, setup$start_var_theta,
                    setup$start_adapt_iter, setup$start_tau_theta)
  cov_ls2_cand = lapply(1:nexp, function(i) {
    AMcov_pool(ntemps, setup$ns2[[i]], setup$start_var_ls2,
               setup$start_adapt_iter, setup$start_tau_ls2)
  })
  S_clust = array(0, c(ntemps, K, p, p))
  for (t in 1:ntemps) for (k in 1:K) S_clust[t, k, , ] = diag(p) * setup$start_var_theta

  count = rep(0, ntemps)
  count_s2 = matrix(0, nexp, ntemps)
  count_temper = matrix(0, ntemps, ntemps)
  llik_trace = rep(0, nmcmc)

  pb <- progress::progress_bar$new(
    format = "  MCMC [:bar] :current/:total (:percent) in :eta [:tick_rate it/sec]",
    total = nmcmc,
    clear = FALSE,
    width = 60
  )

  for (m in 2:nmcmc) {
    m_cur = m
    theta[m, , , ] = theta[m - 1, , , ]
    theta0[m, , ] = theta0[m - 1, , ]
    Sigma0[m, , , ] = Sigma0[m - 1, , , ]
    eta[m, ] = eta[m - 1, ]
    for (i in 1:nexp) {
      delta[[i]][m, , ] = delta[[i]][m - 1, , ]
      log_s2[[i]][m, , ] = log_s2[[i]][m - 1, , ]
      setup$models[[i]] = step_m(setup$models[[i]])
      if (setup$models[[i]]$stochastic) {
        pred_delta[[i]] = eval_clusters(i, th_now(m))
        llik_delta[[i]] = llik_all(i, pred_delta[[i]])
      }
    }

    # ---- memberships: one theta_j at a time (Neal 2000, algorithm 8) ----
    cnt = member_count(m)
    for (i in 1:nexp) {
      for (j in 1:ntheta_vec[i]) {
        for (t in 1:ntemps) {
          cur = delta[[i]][m, t, j]
          cnt[t, cur] = cnt[t, cur] - 1L
          empty = cnt[t, ] == 0
          w = cnt[t, ] + empty * eta[m, t] / (sum(empty) + 1e-9)
          lp = log(w) + itl[t] * llik_delta[[i]][t, , j]
          lp = lp - max(lp)
          new = sample.int(K, 1, prob = exp(lp))
          delta[[i]][m, t, j] = new
          cnt[t, new] = cnt[t, new] + 1L
        }
      }
    }
    ext = extant(m)
    n_ext = rowSums(ext)

    pred_curr = vector("list", nexp)
    llik_curr = vector("list", nexp)
    for (i in 1:nexp) {
      theta_hist[[i]][m, , , ] = hist_from(th_now(m), dl_now(i, m))
      pred_curr[[i]] = pred_from_delta(i, pred_delta[[i]], dl_now(i, m))
    }

    # ---- error variances ----
    for (i in 1:nexp) {
      s2mode = setup$models[[i]]$s2
      if (s2mode == "gibbs") {
        dev_sq = weighted_dev_sq(pred_curr[[i]], setup$ys[[i]], wt[[i]], s2_ind_mat[[i]])
        for (t in 1:ntemps) {
          log_s2[[i]][m, t, ] = gibbs_log_s2(dev_sq[t, ], ny_w[[i]], itl[t],
                                             setup$ig_a[[i]], setup$ig_b[[i]],
                                             setup$sd_lower[[i]], setup$sd_upper[[i]])
          cov_curr[[i]][[t]] = lapply(1:ntheta_vec[i], function(j) cov_ij(i, log_s2[[i]][m, t, ], j))
        }
      } else if (s2mode != "fix") {
        llik_now = llik_given_delta(i, llik_delta[[i]], dl_now(i, m))
        cov_ls2_cand[[i]] = update_m(cov_ls2_cand[[i]], log_s2[[i]], m)
        ls2_cand = gen_cand(cov_ls2_cand[[i]], log_s2[[i]], m)
        ls2_curr = matrix(log_s2[[i]][m - 1, , ], ntemps, setup$ns2[[i]])
        llik_candi = matrix(0, ntemps, ntheta_vec[i])
        cov_candi = vector("list", ntemps)
        for (t in 1:ntemps) {
          ls2_t = pmin(ls2_cand[t, ], log(1e100))
          cov_candi[[t]] = lapply(1:ntheta_vec[i], function(j) cov_ij(i, ls2_t, j))
          for (j in 1:ntheta_vec[i]) {
            llik_candi[t, j] = llik_ij(i, pred_curr[[i]][t, ], cov_candi[[t]][[j]], j)
          }
        }
        alpha_s2 = itl * (rowSums(llik_candi) - rowSums(llik_now))
        alpha_s2 = alpha_s2 + itl * s2_kern_sum(setup$s2_prior_kern[[i]], ls2_cand, setup$ig_a[[i]], setup$ig_b[[i]])
        alpha_s2 = alpha_s2 + itl * ls2_rowsum(ls2_cand)
        alpha_s2 = alpha_s2 - itl * s2_kern_sum(setup$s2_prior_kern[[i]], ls2_curr, setup$ig_a[[i]], setup$ig_b[[i]])
        alpha_s2 = alpha_s2 - itl * ls2_rowsum(ls2_curr)
        alpha_s2[!ls2_in_bounds(ls2_cand, setup$sd_lower[[i]], setup$sd_upper[[i]])] = -Inf
        for (t in which(log(stats::runif(ntemps)) < alpha_s2)) {
          count_s2[i, t] = count_s2[i, t] + 1
          log_s2[[i]][m, t, ] = ls2_cand[t, ]
          cov_curr[[i]][[t]] = cov_candi[[t]]
          cov_ls2_cand[[i]]$count_100[t] = cov_ls2_cand[[i]]$count_100[t] + 1
        }
        cov_ls2_cand[[i]] = update_tau(cov_ls2_cand[[i]], m)
      }
    }

    # ---- Gibbs update theta0 from the occupied clusters ----
    th = th_now(m)
    for (t in 1:ntemps) {
      occ = matrix(th[t, ext[t, ], ], ncol = p)
      prec = n_ext[t] * itl[t] * Sigma0_inv[t, , ] + theta0_prior_prec
      cc = chol2inv(chol(prec))
      dd = itl[t] * Sigma0_inv[t, , ] %*% colSums(occ) + theta0_prior_prec %*% theta0_prior_mean
      theta0[m, t, ] = rmvn_constrained(t(cc %*% dd), t(chol(cc)), setup)
    }

    # ---- Gibbs update Sigma0 ----
    for (t in 1:ntemps) {
      occ = matrix(th[t, ext[t, ], ], ncol = p)
      d = sweep(occ, 2, theta0[m, t, ])
      S_t = rinvwishart(Sigma0_prior_df + n_ext[t] * itl[t], Sigma0_prior_scale + itl[t] * crossprod(d))
      Sigma0[m, t, , ] = S_t
      R = chol(S_t)
      Sigma0_inv[t, , ] = chol2inv(R)
      Sigma0_ldet[t] = 2 * sum(log(diag(R)))
    }

    # ---- cluster values: adaptive Metropolis on occupied, prior draw on empty ----
    am = am_hist_update(am, theta_hist, m)
    if (isTRUE(am$adapted)) {
      for (t in 1:ntemps) for (k in 1:K) {
        members = list()
        for (i in 1:nexp) {
          for (j in which(delta[[i]][m, t, ] == k)) {
            members[[length(members) + 1]] = list(mu = am$mu[[i]][t, j, ], cov = am$cov[[i]][t, j, , ])
          }
        }
        if (length(members) == 0) next
        mus = do.call(rbind, lapply(members, `[[`, "mu"))
        mu_k = colMeans(mus)
        cov_k = Reduce(`+`, lapply(members, `[[`, "cov")) / length(members) +
          crossprod(sweep(mus, 2, mu_k)) / length(members)
        S_clust[t, k, , ] = am$AM_SCALAR * exp(am$tau[[1]][t, 1]) * (cov_k + diag(p) * am$eps)
      }
    }
    cnt = member_count(m)

    cand = th
    for (t in 1:ntemps) for (k in 1:K) {
      if (ext[t, k]) {
        S = S_clust[t, k, , ]
        if (setup$scale_proposal) S = S / max(cnt[t, k], 1)
        cand[t, k, ] = th[t, k, ] + t(chol(S)) %*% stats::rnorm(p)
      }
    }
    good = matrix(theta_ok(matrix(cand, ntemps * K, p), setup), ntemps, K) & ext
    # empty clusters get a fresh draw from the hierarchical distribution
    for (t in 1:ntemps) {
      if (any(!ext[t, ])) {
        fresh = rmvn_constrained(matrix(theta0[m, t, ], sum(!ext[t, ]), p, byrow = TRUE),
                                 t(chol(Sigma0[m, t, , ])), setup)
        cand[t, !ext[t, ], ] = fresh
      }
    }
    eval_arr = th
    for (t in 1:ntemps) for (k in 1:K) {
      if (good[t, k] || !ext[t, k]) eval_arr[t, k, ] = cand[t, k, ]
    }

    pred_delta_cand = vector("list", nexp)
    llik_clust_cand = matrix(0, ntemps, K)
    llik_clust_curr = matrix(0, ntemps, K)
    for (i in 1:nexp) {
      pred_delta_cand[[i]] = eval_clusters(i, eval_arr)
      dl = dl_now(i, m)
      pc = pred_from_delta(i, pred_delta_cand[[i]], dl)
      for (t in 1:ntemps) for (j in 1:ntheta_vec[i]) {
        k = dl[t, j]
        llik_clust_cand[t, k] = llik_clust_cand[t, k] + llik_ij(i, pc[t, ], cov_curr[[i]][[t]][[j]], j)
        llik_clust_curr[t, k] = llik_clust_curr[t, k] + llik_ij(i, pred_curr[[i]][t, ], cov_curr[[i]][[t]][[j]], j)
      }
    }
    alpha = matrix(-Inf, ntemps, K)
    lp_diff = lprior_clusters(eval_arr) - lprior_clusters(th)
    alpha[good] = (itl * (llik_clust_cand - llik_clust_curr + lp_diff))[good]
    accept = (log(matrix(stats::runif(ntemps * K), ntemps)) < alpha) & good
    update = accept | !ext
    for (t in 1:ntemps) for (k in which(update[t, ])) theta[m, t, k, ] = eval_arr[t, k, ]
    for (i in 1:nexp) {
      for (t in 1:ntemps) for (k in which(update[t, ])) {
        pred_delta[[i]][t, k, ] = pred_delta_cand[[i]][t, k, ]
      }
      theta_hist[[i]][m, , , ] = hist_from(th_now(m), dl_now(i, m))
      llik_delta[[i]] = llik_all(i, pred_delta[[i]])
      pred_curr[[i]] = pred_from_delta(i, pred_delta[[i]], dl_now(i, m))
      llik_curr[[i]] = llik_given_delta(i, llik_delta[[i]], dl_now(i, m))
    }
    count = count + rowSums(accept)

    # ---- concentration eta (Escobar & West 1995) ----
    for (t in 1:ntemps) {
      g = stats::rbeta(1, eta[m, t] + 1, ndat)
      aa = eta_shape + n_ext[t]
      bb = eta_rate - log(g)
      pi_hi = (aa - 1) / (ndat * bb + aa - 1)
      shape = if (stats::runif(1) < pi_hi) aa else aa - 1
      eta[m, t] = stats::rgamma(1, shape = shape, rate = bb)
    }

    # ---- tempering swaps ----
    if ((m > setup$start_temper) && (ntemps > 1)) {
      for (sweep_k in 1:setup$nswap) {
        ext = extant(m)
        lpc = lprior_clusters(th_now(m))
        tempered = rowSums(lpc * ext)
        for (i in 1:nexp) {
          ls2 = matrix(log_s2[[i]][m, , ], ntemps, setup$ns2[[i]])
          tempered = tempered + rowSums(llik_curr[[i]]) +
            s2_kern_sum(setup$s2_prior_kern[[i]], ls2, setup$ig_a[[i]], setup$ig_b[[i]])
        }
        sw = matrix(sample(1:ntemps, 2 * setup$nswap_per, FALSE), 2)
        sw_alpha = (itl[sw[2, ]] - itl[sw[1, ]]) * (tempered[sw[1, ]] - tempered[sw[2, ]])
        for (s in which(log(stats::runif(setup$nswap_per)) < sw_alpha)) {
          a = sw[1, s]
          b = sw[2, s]
          count_temper[a, b] = count_temper[a, b] + 1
          theta[m, c(a, b), , ] = theta[m, c(b, a), , ]
          theta0[m, c(a, b), ] = theta0[m, c(b, a), ]
          Sigma0[m, c(a, b), , ] = Sigma0[m, c(b, a), , ]
          Sigma0_inv = swap_rows(Sigma0_inv, a, b)
          Sigma0_ldet = swap_rows(Sigma0_ldet, a, b)
          eta[m, c(a, b)] = eta[m, c(b, a)]
          for (i in 1:nexp) {
            delta[[i]][m, c(a, b), ] = delta[[i]][m, c(b, a), ]
            theta_hist[[i]][m, c(a, b), , ] = theta_hist[[i]][m, c(b, a), , ]
            log_s2[[i]][m, c(a, b), ] = log_s2[[i]][m, c(b, a), ]
            pred_delta[[i]] = swap_rows(pred_delta[[i]], a, b)
            llik_delta[[i]] = swap_rows(llik_delta[[i]], a, b)
            pred_curr[[i]] = swap_rows(pred_curr[[i]], a, b)
            llik_curr[[i]] = swap_rows(llik_curr[[i]], a, b)
            cov_curr[[i]][c(a, b)] = cov_curr[[i]][c(b, a)]
          }
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
    lapply(1:ntheta_vec[i], function(j) to_native(theta_hist[[i]][, 1, j, ]))
  })
  parent = matrix(0, nmcmc, p)
  for (mm in 1:nmcmc) {
    parent[mm, ] = rmvn_constrained(theta0[mm, 1, ], t(chol(Sigma0[mm, 1, , ])), setup)
  }
  count_temper = count_temper + t(count_temper) - diag(diag(count_temper), ntemps)

  list(
    theta = theta,
    theta_hist = theta_hist,
    delta = delta,
    eta = eta,
    theta0 = theta0,
    Sigma0 = Sigma0,
    s2 = s2,
    count = count,
    count_s2 = count_s2,
    count_temper = count_temper,
    pred_curr = pred_curr,
    llik = llik_trace,
    nclustmax = K,
    theta_native = theta_native,
    theta0_native = to_native(theta0[, 1, ]),
    theta_parent_native = to_native(parent)
  )
}
