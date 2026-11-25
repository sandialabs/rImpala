# Internal helpers shared by calibPool(), calibHier() and calibClust().


# Marginal likelihood covariance for the mvBayes family of emulators. `inds`
# picks out the observations `s2vec` belongs to, so hierarchical and clustered
# calibration can build the covariance of one theta's subset of yobs.
mv_lik_cov_inv <- function(obj, s2vec, inds = NULL) {
  if (is.null(inds)) {
    inds = seq_along(s2vec)
  }
  basis = obj$basis[inds, , drop = FALSE]
  Sigma = cor2cov(obj$meas_error_corr[inds, inds, drop = FALSE], sqrt(s2vec))
  mat = Sigma + obj$trunc_error_var[inds, inds, drop = FALSE] +
    obj$discrep_cov[inds, inds, drop = FALSE] +
    basis %*% diag(obj$emu_vars, nrow = obj$npc) %*% t(basis)
  R <- tryCatch({
    chol(mat)
  }, error = function(e) {
    chol(mat + 0.0001 * diag(nrow(mat)))
  })
  ldet = 2 * sum(log(diag(R)))
  inv = chol2inv(R)
  list(inv = inv, ldet = ldet)
}


# W^(1/2) %*% inv %*% W^(1/2) for a vector of observation weights.
weight_inv <- function(inv, wt) {
  if (is.null(wt)) {
    return(inv)
  }
  sw = sqrt(wt)
  inv * outer(sw, sw)
}


# Observation weights for experiment i; setups built before weights existed
# carry none, which means every weight is 1.
get_wt <- function(setup, i) {
  wt = setup$wt[[i]]
  if (is.null(wt)) {
    wt = rep(1, length(setup$ys[[i]]))
  }
  wt
}


# Non-unit weights only mean something to a model whose llik/lik_cov_inv apply
# them. Refuse rather than silently fit an unweighted likelihood.
check_wt_support <- function(setup) {
  for (i in seq_len(setup$nexp)) {
    wt = get_wt(setup, i)
    if (any(wt != 1) && !isTRUE(setup$models[[i]]$supports_wt)) {
      cli::cli_abort(c(
        "Experiment {i} has observation weights, but its model ({.cls {class(setup$models[[i]])[1]}}) does not apply them.",
        i = "Weights are implemented for {.fn ModelBassPca_func}.",
        i = "A custom model can opt in by setting {.code supports_wt = TRUE} and using the {.arg wt} argument of its {.fn llik} and {.fn lik_cov_inv} methods."
      ))
    }
  }
  invisible(NULL)
}


# Gibbs draw of the log error variances for one temperature. The likelihood is
# a product of N(r_k; 0, s2)^wt_k over observations, which stays conjugate to
# the inverse gamma prior for any weights, so the weighted residual sum of
# squares and weighted observation count simply replace the unweighted ones.
# With unit weights this is exactly the unweighted update.
gibbs_log_s2 <- function(dev_w, ny_w, itl, ig_a, ig_b, sd_lower, sd_upper) {
  shape = itl * (ny_w / 2 + ig_a + 1) - 1
  scale = 1 / (itl * (ig_b + dev_w / 2))
  log(rig_bounded(shape, scale, sd_lower, sd_upper))
}


# Weighted squared residuals summed per variance group: (ntemps x ns2).
# `pred` is (ntemps x ny), `s2_ind_mat` the (ny x ns2) indicator.
weighted_dev_sq <- function(pred, yobs, wt, s2_ind_mat) {
  pred = matrix(pred, ncol = length(yobs))
  dev = sweep(pred, 2, yobs)^2
  sweep(dev, 2, wt, "*") %*% s2_ind_mat
}


# Call the setup's constraint function on a named list of native-scale
# parameters. Functions that accept a third argument (or `...`) also get the
# model constants, so a constraint can involve parameters held constant (see
# constraints_ptw); two-argument functions are called as before.
check_constraints <- function(setup, x) {
  f = setup$checkConstraints
  fa = names(formals(f))
  if (length(fa) >= 3 || "..." %in% fa) {
    f(x, setup$bounds, setup$constants)
  } else {
    f(x, setup$bounds)
  }
}


# TRUE/FALSE per row of an (n x p) 0-1 scale matrix: satisfies the constraints.
theta_ok <- function(x, setup) {
  x = matrix(x, ncol = setup$p)
  as.logical(check_constraints(setup, tran_unif(x, setup$bounds_mat, names(setup$bounds))))
}


# Uniform draws on the unit cube that satisfy the constraints, (n x p).
sample_unif_constrained <- function(n, setup, maxiter = 10000) {
  x = matrix(stats::runif(n * setup$p), n, setup$p)
  good = theta_ok(x, setup)
  it = 0
  while (any(!good)) {
    it = it + 1
    if (it > maxiter) {
      cli::cli_abort("Could not find starting values satisfying the constraints after {maxiter} tries.")
    }
    nbad = sum(!good)
    x[!good, ] = matrix(stats::runif(nbad * setup$p), nbad, setup$p)
    good[!good] = theta_ok(x[!good, , drop = FALSE], setup)
  }
  x
}


# One multivariate normal draw per row of `means` (n x p), all sharing the
# covariance whose lower Cholesky factor is `L`, resampled until each row
# satisfies the constraints (a truncated normal).
rmvn_constrained <- function(means, L, setup, maxiter = 10000) {
  means = matrix(means, ncol = setup$p)
  n = nrow(means)
  draw = function(k) {
    means[k, , drop = FALSE] +
      t(L %*% matrix(stats::rnorm(setup$p * length(k)), setup$p, length(k)))
  }
  x = draw(seq_len(n))
  good = theta_ok(x, setup)
  it = 0
  while (any(!good)) {
    it = it + 1
    if (it > maxiter) {
      cli::cli_abort(c(
        "Could not draw parameters satisfying the constraints after {maxiter} tries.",
        i = "The hierarchical distribution N(theta0, Sigma0) puts almost no mass inside the constraints."
      ))
    }
    bad = which(!good)
    x[bad, ] = draw(bad)
    good[bad] = theta_ok(x[bad, , drop = FALSE], setup)
  }
  x
}


# Inverse-Wishart draw with the same (df, scale) convention as scipy and the
# python sampler: E[Sigma] = scale / (df - p - 1).
rinvwishart <- function(df, scale) {
  W = stats::rWishart(1, df, chol2inv(chol(scale)))[, , 1]
  chol2inv(chol(W))
}


# Log density of an inverse Wishart, for the tempering swap.
dinvwishart_log <- function(Sigma, df, scale) {
  p = nrow(Sigma)
  ld_scale = 2 * sum(log(diag(chol(scale))))
  ld_Sigma = 2 * sum(log(diag(chol(Sigma))))
  0.5 * df * ld_scale - 0.5 * df * p * log(2) - lmvgamma(df / 2, p) -
    0.5 * (df + p + 1) * ld_Sigma - 0.5 * sum(diag(scale %*% chol2inv(chol(Sigma))))
}


lmvgamma <- function(a, p) {
  p * (p - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(p)) / 2))
}


# Log N(x; mean, Sigma) for each row of x (n x p), given Sigma's inverse and
# log determinant.
dmvn_rows <- function(x, mean, prec, ldet) {
  x = matrix(x, ncol = length(mean))
  d = sweep(x, 2, mean)
  -0.5 * ncol(x) * log(2 * pi) - 0.5 * ldet - 0.5 * rowSums((d %*% prec) * d)
}


# Evaluate an emulator at every (temperature, theta) row and keep, for each
# observation, the prediction made with the theta that observation belongs to.
# `theta` is (ntemps x ntheta x p) on the 0-1 scale; `theta_ind` maps the ny
# observations to 1..ntheta. Returns (ntemps x ny).
eval_by_theta <- function(model, theta, theta_ind, setup) {
  ntemps = dim(theta)[1]
  ntheta = dim(theta)[2]
  # row r of `mat` is temperature ((r - 1) %% ntemps) + 1 and theta ((r - 1) %/% ntemps) + 1
  mat = matrix(theta, ntemps * ntheta, setup$p)
  full = evalm(model, tran_unif(mat, setup$bounds_mat, names(setup$bounds)), TRUE)
  ny = length(theta_ind)
  full = matrix(full, ntemps * ntheta, ny)
  rows = rep(seq_len(ntemps), times = ny) + rep((theta_ind - 1) * ntemps, each = ntemps)
  cols = rep(seq_len(ny), each = ntemps)
  matrix(full[cbind(rows, cols)], ntemps, ny)
}


# Checks shared by calibHier() and calibClust().
check_hier_setup <- function(setup, fname) {
  if (setup$nexp == 0) {
    cli::cli_abort("Add experiments with {.fn addVecExperiments} before calling {.fn {fname}}.")
  }
  if (is.null(setup$theta0_prior_mean)) {
    cli::cli_abort("Call {.fn setHierPriors} before {.fn {fname}}.")
  }
  if (length(setup$theta_fixed) > 0) {
    cli::cli_abort(c(
      "{.fn fixTheta} is not supported by {.fn {fname}} yet.",
      i = "Use {.fn calibPool} for conditional inference, or remove the fixed parameters with {.fn unfixTheta}."
    ))
  }
  if (length(setup$theta_prior) > 0) {
    cli::cli_warn(c(
      "Priors added with {.fn addThetaPrior} are ignored by {.fn {fname}}.",
      i = "The parameters get the hierarchical prior set by {.fn setHierPriors} instead."
    ))
  }
  check_wt_support(setup)
  invisible(NULL)
}


# ntemps x p matrix of starting values for theta0.
theta0_start_values <- function(setup) {
  if (!is.null(setup$theta0_start)) {
    st = matrix(setup$theta0_start, ncol = setup$p)
    if (nrow(st) == 1) {
      st = matrix(st, setup$ntemps, setup$p, byrow = TRUE)
    }
    if (nrow(st) != setup$ntemps) {
      cli::cli_abort("{.arg theta0_start} has {nrow(st)} rows but there are {setup$ntemps} temperatures.")
    }
    if (!all(theta_ok(st, setup))) {
      cli::cli_abort("{.arg theta0_start} does not satisfy the constraints.")
    }
    return(st)
  }
  sample_unif_constrained(setup$ntemps, setup)
}


# Running-moment adaptive Metropolis state for a list of (nmcmc x ntemps x
# ntheta_i x p) histories, one proposal covariance per (temperature, theta).
am_hist_init <- function(ntemps, ntheta, p, start_var, start_adapt_iter, tau_start) {
  nexp = length(ntheta)
  S = lapply(ntheta, function(n) {
    a = array(0, c(ntemps, n, p, p))
    for (t in 1:ntemps) for (j in 1:n) a[t, j, , ] = diag(p) * start_var
    a
  })
  list(
    eps = 1e-12,
    AM_SCALAR = 2.4^2 / p,
    tau = lapply(ntheta, function(n) matrix(tau_start, ntemps, n)),
    S = S,
    cov = lapply(ntheta, function(n) array(0, c(ntemps, n, p, p))),
    mu = lapply(ntheta, function(n) array(0, c(ntemps, n, p))),
    nexp = nexp,
    ntheta = ntheta,
    ntemps = ntemps,
    p = p,
    start_adapt_iter = start_adapt_iter,
    count_100 = lapply(ntheta, function(n) matrix(0, ntemps, n))
  )
}


# Update the running mean/covariance with draw m - 1. Returns the updated
# object with `adapted = TRUE` once proposal covariances are being learned.
am_hist_update <- function(obj, x, m) {
  p = obj$p
  if (m > obj$start_adapt_iter) {
    for (i in 1:obj$nexp) {
      xprev = array(x[[i]][m - 1, , , ], c(obj$ntemps, obj$ntheta[i], p))
      obj$mu[[i]] = obj$mu[[i]] + (xprev - obj$mu[[i]]) / m
      d = xprev - obj$mu[[i]]
      for (t in 1:obj$ntemps) for (j in 1:obj$ntheta[i]) {
        obj$cov[[i]][t, j, , ] = ((m - 1) / m) * obj$cov[[i]][t, j, , ] +
          ((m - 1) / (m * m)) * tcrossprod(d[t, j, ])
      }
    }
    obj$adapted = TRUE
  } else if (m == obj$start_adapt_iter) {
    for (i in 1:obj$nexp) {
      xh = array(x[[i]][1:m, , , ], c(m, obj$ntemps, obj$ntheta[i], p))
      for (t in 1:obj$ntemps) for (j in 1:obj$ntheta[i]) {
        xm = matrix(xh[, t, j, ], m, p)
        obj$mu[[i]][t, j, ] = colMeans(xm)
        obj$cov[[i]][t, j, , ] = stats::cov(xm)
      }
    }
    obj$adapted = TRUE
  }
  obj
}


# Proposal covariances from the running covariances (hierarchical sampler).
am_hist_scale <- function(obj) {
  if (!isTRUE(obj$adapted)) {
    return(obj)
  }
  for (i in 1:obj$nexp) for (t in 1:obj$ntemps) for (j in 1:obj$ntheta[i]) {
    obj$S[[i]][t, j, , ] = obj$AM_SCALAR * exp(obj$tau[[i]][t, j]) *
      (obj$cov[[i]][t, j, , ] + diag(obj$p) * obj$eps)
  }
  obj
}


am_hist_tau <- function(obj, m) {
  if ((m %% 100 == 0) && (m > obj$start_adapt_iter)) {
    delta = min(0.5, 5 / sqrt(m + 1))
    for (i in 1:obj$nexp) {
      low = obj$count_100[[i]] < 23
      high = obj$count_100[[i]] > 23
      obj$tau[[i]][low] = obj$tau[[i]][low] - delta
      obj$tau[[i]][high] = obj$tau[[i]][high] + delta
      obj$count_100[[i]][] = 0
    }
  }
  obj
}


# Swap two temperature slices of an array along its first (non-mcmc) axis.
swap_rows <- function(x, a, b) {
  if (is.null(dim(x))) {
    x[c(a, b)] = x[c(b, a)]
  } else if (length(dim(x)) == 2) {
    x[c(a, b), ] = x[c(b, a), ]
  } else if (length(dim(x)) == 3) {
    x[c(a, b), , ] = x[c(b, a), , ]
  } else {
    x[c(a, b), , , ] = x[c(b, a), , , ]
  }
  x
}
