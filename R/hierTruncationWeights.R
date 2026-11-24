#' @title Diagnose the bounds-truncation contraction of a hierarchical fit
#'
#' @description [calibHier()] is a correct sampler for the joint density
#'   \deqn{p(\theta_0) p(\Sigma_0) \prod_j N(\theta_j; \theta_0, \Sigma_0)
#'     1\{\theta_j \in B\} \times \mathrm{likelihood},}
#'   but integrating out \eqn{\theta} leaves the effective hyperprior equal to
#'   the stated one times \eqn{\prod_j P(\theta_j \in B \mid \theta_0,
#'   \Sigma_0)}, because the truncation to the bounds \eqn{B = [0, 1]^p} is
#'   applied without its normalising constant. That factor rewards a central
#'   \eqn{\theta_0} and a small \eqn{\Sigma_0}, so the fitted hierarchy looks
#'   tighter than the prior that was written down. The intended model --
#'   \eqn{\theta_j} truncated normal with its normalising constant included --
#'   is recovered by reweighting the saved cold-chain draws with
#'   \deqn{w \propto P(B \mid \theta_0, \Sigma_0)^{-n},}
#'   where \eqn{n} is `sum(setup$ntheta)`, the total number of
#'   experiment-level \eqn{\theta_j}. This function returns
#'   those weights together with the raw and reweighted summaries of
#'   \eqn{\theta_0} and \eqn{\Sigma_0}, so a finished run can be checked (and,
#'   where feasible, corrected) after the fact.
#'
#' @param setup the `CalibSetup` object passed to [calibHier()]
#' @param out the list [calibHier()] returned
#' @param mcmc_use iterations to use, e.g. `seq(nburn + 1, nmcmc)`
#'   (default: `NULL`, the second half of the chain)
#' @param temperature temperature index, 1-based; `1` (the default) is the cold
#'   chain
#'
#' @return A list with
#' - `weights`: normalised weights, one per used draw
#' - `ess`, `ess_frac`: effective sample size \eqn{1 / \sum w^2} of the
#'   reweighting, and that divided by the number of draws used
#' - `log_box_prob`: \eqn{\log P(B \mid \theta_0, \Sigma_0)} at each used draw
#' - `n_theta`: \eqn{n}, the exponent of the correction
#' - `tilt_nats`: \eqn{n \,|\mathrm{median} \log P(B)|}, the size of the
#'   correction in nats
#' - `theta0_mean_raw`, `theta0_mean_reweighted`: length `p` means of
#'   \eqn{\theta_0} (0-1 scale)
#' - `theta0_sd_raw`, `theta0_sd_reweighted`: length `p` standard deviations of
#'   \eqn{\theta_0}, both with denominator \eqn{n} rather than \eqn{n - 1}, so
#'   that the two are directly comparable
#' - `sigma0_diag_median_raw`, `sigma0_diag_median_reweighted`: length `p`
#'   medians of the diagonal of \eqn{\Sigma_0}, the reweighted one a
#'   [weightedQuantile()]
#'
#' @details Read `ess_frac` first:
#'
#' - near 1: the bounds are not binding for this run and the uncorrected output
#'   needs no adjustment;
#' - moderate: the reweighted summaries are usable;
#' - small: the output is materially contracted *and* reweighting cannot repair
#'   it, because the chain never visited the region the corrected posterior
#'   wants. Treat the reported hyperprior as an approximation, or reduce the
#'   number of parameters sitting against a bound.
#'
#'   `tilt_nats` grows with both the number of experiments and the dimension,
#'   and is only large when parameters sit near the bounds. The contraction is a
#'   weak-data effect: with every \eqn{\theta_j} identified by its own data the
#'   parent spread is recovered well and only the coordinates nearest a bound
#'   are pulled inward.
#'
#'   \eqn{P(B \mid \theta_0, \Sigma_0)} is computed in closed form, as a product
#'   of normal CDFs, when \eqn{\Sigma_0} is diagonal, and with
#'   `mvtnorm::pmvnorm()` otherwise. The weights account for the bounds only:
#'   with a `constraint_func` other than [cf_bounds()] the true normalising
#'   constant is smaller than the box probability, so the correction is partial
#'   and a warning is issued rather than silently under-correcting.
#'
#' @seealso [calibHier()], [setHierPriors()], [weightedQuantile()]
#'
#' @export
#'
hierTruncationWeights <- function(setup, out, mcmc_use = NULL,
                                  temperature = 1) {
  if (is.null(setup$theta0_prior_mean)) {
    cli::cli_abort(c(
      "{.arg setup} does not carry hierarchical hyperpriors.",
      i = "Pass the {.cls CalibSetup} that {.fn calibHier} was run on, after {.fn setHierPriors}."
    ))
  }
  if (is.null(out$theta0) || is.null(out$Sigma0)) {
    cli::cli_abort("{.arg out} must have {.field theta0} and {.field Sigma0}, as returned by {.fn calibHier}.")
  }

  theta0 = out$theta0
  Sigma0 = out$Sigma0
  nmcmc = dim(theta0)[1]
  ntemps = dim(theta0)[2]
  p = dim(theta0)[3]

  if (length(temperature) != 1 || !is.finite(temperature) ||
      temperature < 1 || temperature > ntemps ||
      temperature != as.integer(temperature)) {
    cli::cli_abort("{.arg temperature} must be a whole number between 1 and {ntemps}; 1 is the cold chain.")
  }
  if (is.null(mcmc_use)) {
    mcmc_use = seq.int(floor(nmcmc / 2) + 1, nmcmc)
  }
  mcmc_use = as.integer(mcmc_use)
  if (length(mcmc_use) < 1 || anyNA(mcmc_use) || any(mcmc_use < 1) ||
      any(mcmc_use > nmcmc)) {
    cli::cli_abort("{.arg mcmc_use} must be a non-empty set of iteration indices between 1 and {nmcmc}.")
  }

  n_theta = sum(setup$ntheta)
  if (length(n_theta) != 1 || !is.finite(n_theta) || n_theta < 1) {
    cli::cli_abort("{.arg setup} has no experiment-level {.field theta}; was {.fn addVecExperiments} called with {.arg theta_ind}?")
  }

  if (!identical(setup$checkConstraints, cf_bounds)) {
    cli::cli_warn(c(
      "{.fn hierTruncationWeights} accounts for the bounds only.",
      i = "This setup uses a custom {.arg constraint_func}, so the correction is partial."
    ))
  }

  th0 = matrix(theta0[mcmc_use, temperature, ], ncol = p)
  nd = length(mcmc_use)
  S0 = lapply(mcmc_use, function(m) matrix(Sigma0[m, temperature, , ], p, p))

  lp = rep(NA_real_, nd)
  for (j in seq_len(nd)) {
    lp[j] = tryCatch(log_box_prob(th0[j, ], S0[[j]]),
                     error = function(e) NA_real_)
  }
  ok = is.finite(lp)
  if (!all(ok)) {
    cli::cli_warn(
      "{sum(!ok)} of {nd} draw{?s} gave an unusable box probability and {?was/were} dropped from the diagnostic."
    )
  }
  if (!any(ok)) {
    cli::cli_abort("No usable draws: every box probability failed. Are the {.field Sigma0} draws positive definite?")
  }
  th0 = matrix(th0[ok, ], ncol = p)
  S0 = S0[ok]
  lp = lp[ok]

  lw = -n_theta * lp
  w = exp(lw - max(lw))
  w = w / sum(w)
  ess = 1 / sum(w^2)

  flat = rep(1 / length(w), length(w))
  # p x nd then transposed, so that p = 1 also gives an nd x 1 matrix
  s0diag = t(matrix(vapply(S0, diag, numeric(p)), nrow = p))
  dimnames(s0diag) = NULL

  mean_raw = as.numeric(crossprod(flat, th0))
  mean_rw = as.numeric(crossprod(w, th0))

  list(
    weights = w,
    ess = as.numeric(ess),
    ess_frac = as.numeric(ess / length(w)),
    log_box_prob = lp,
    n_theta = as.integer(n_theta),
    tilt_nats = as.numeric(n_theta * abs(stats::median(lp))),
    theta0_mean_raw = mean_raw,
    theta0_mean_reweighted = mean_rw,
    theta0_sd_raw = weighted_sd(th0, flat, mean_raw),
    theta0_sd_reweighted = weighted_sd(th0, w, mean_rw),
    sigma0_diag_median_raw = apply(s0diag, 2, stats::median),
    sigma0_diag_median_reweighted = vapply(
      seq_len(p),
      function(k) weightedQuantile(s0diag[, k], 0.5, w),
      numeric(1)
    )
  )
}


#' @title Weighted quantiles
#'
#' @description Quantiles of a one-dimensional sample carrying unequal weights,
#'   as produced by [hierTruncationWeights()]. Uses the weighted analogue of
#'   `type = 5` plotting positions: the quantile level of the `i`-th ordered
#'   value is its cumulative weight less half its own weight, and levels in
#'   between are interpolated linearly, clamped to the sample range.
#'
#' @param x numeric vector of values
#' @param probs numeric vector of probabilities in `[0, 1]`
#' @param weights numeric vector of non-negative weights, the same length as
#'   `x`; need not sum to 1 (default: `NULL`, equal weights)
#'
#' @return A numeric vector of quantiles, one per entry of `probs`.
#'
#' @seealso [hierTruncationWeights()]
#'
#' @export
#'
weightedQuantile <- function(x, probs, weights = NULL) {
  x = as.numeric(x)
  if (is.null(weights)) weights = rep(1, length(x))
  weights = as.numeric(weights)
  if (length(weights) != length(x)) {
    cli::cli_abort("{.arg weights} must be the same length as {.arg x}.")
  }
  keep = is.finite(x) & is.finite(weights) & weights > 0
  if (!any(keep)) {
    cli::cli_abort("{.arg x} and {.arg weights} have no finite entry with positive weight.")
  }
  x = x[keep]
  weights = weights[keep]
  if (any(probs < 0 | probs > 1)) {
    cli::cli_abort("{.arg probs} must lie in {.code [0, 1]}.")
  }
  if (length(x) == 1) return(rep(x, length(probs)))

  ord = order(x)
  x = x[ord]
  weights = weights[ord]
  cw = cumsum(weights) - 0.5 * weights
  cw = cw / sum(weights)
  stats::approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}


# log P(theta in [0, 1]^p) for theta ~ N(mean, cov), the normalising constant
# calibHier drops. Closed form (a product of normal CDFs) when cov is diagonal,
# mvtnorm::pmvnorm (Genz's algorithm) otherwise. pmvnorm draws from the R RNG,
# so the stream is restored on exit to keep this diagnostic side-effect free.
log_box_prob <- function(mean, cov) {
  mean = as.numeric(mean)
  cov = matrix(cov, length(mean), length(mean))
  p = length(mean)
  floor_log = log(1e-300)

  if (p == 1) {
    sd = sqrt(cov[1, 1])
    val = stats::pnorm((1 - mean) / sd) - stats::pnorm(-mean / sd)
    return(max(log(max(val, 0)), floor_log))
  }

  off = cov - diag(diag(cov), p)
  if (max(abs(off)) <= 1e-12 * max(max(abs(cov)), 1)) {
    sd = sqrt(diag(cov))
    val = stats::pnorm((1 - mean) / sd) - stats::pnorm(-mean / sd)
    return(sum(log(pmax(val, 1e-300))))
  }

  seed = if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  on.exit({
    if (is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", seed, envir = globalenv())
    }
  }, add = TRUE)
  val = mvtnorm::pmvnorm(
    lower = rep(0, p), upper = rep(1, p), mean = mean, sigma = cov,
    algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 1e-8, releps = 0)
  )
  max(log(max(as.numeric(val), 0)), floor_log)
}


# sd with denominator n (weights summing to 1), so that the uniformly weighted
# "raw" summary and the reweighted one are computed the same way.
weighted_sd <- function(x, w, mean) {
  sqrt(as.numeric(crossprod(w, sweep(x, 2, mean)^2)))
}
