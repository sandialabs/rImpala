# The bounds-truncation diagnostic for calibHier. calibHier applies the
# truncation of theta_j ~ N(theta0, Sigma0) to the bounds without its
# normalising constant, so the effective hyperprior is tilted by
# prod_j P(B | theta0, Sigma0); hierTruncationWeights() reweights the saved
# draws by that factor to the power -n. No MCMC is needed to check the weights
# themselves: only setup$ntheta and the saved theta0 / Sigma0 are read.

trunc_setup <- function(p = 2, ng = 3, ny_per = 4, constraint_func = cf_bounds) {
  bounds <- list()
  for (k in seq_len(p)) bounds[[paste0("t_", k)]] <- c(0, 1)
  setup <- CalibSetup(bounds, constraint_func)
  ny <- ng * ny_per
  group <- rep(seq_len(ng), each = ny_per)
  A <- matrix(1, ny, p)
  setup <- addVecExperiments(setup, rep(0, ny), StubModel(A),
                             sd_est = rep(0.1, ng), s2_df = rep(2, ng),
                             s2_ind = group, theta_ind = group)
  setHierPriors(setup)
}

# A calibHier-shaped output holding the given cold-chain draws, with a second
# (hot) temperature that must be ignored.
trunc_out <- function(th0, s0diag, ntemps = 2) {
  nd <- nrow(th0)
  p <- ncol(th0)
  theta0 <- array(0, c(nd, ntemps, p))
  theta0[, 1, ] <- th0
  Sigma0 <- array(0, c(nd, ntemps, p, p))
  for (j in seq_len(nd)) {
    Sigma0[j, 1, , ] <- diag(s0diag[j, ], p)
    if (ntemps > 1) Sigma0[j, 2, , ] <- diag(p)
  }
  list(theta0 = theta0, Sigma0 = Sigma0)
}


test_that("weights match the closed form for a diagonal Sigma0", {
  setup <- trunc_setup(p = 2, ng = 3)
  n <- sum(setup$ntheta)
  expect_equal(n, 3)
  set.seed(0)
  nd <- 50
  th0 <- matrix(stats::runif(nd * 2, 0.1, 0.9), nd, 2)
  s0 <- matrix(stats::runif(nd * 2, 0.02, 0.5), nd, 2)^2
  res <- hierTruncationWeights(setup, trunc_out(th0, s0), mcmc_use = 1:nd)

  sd <- sqrt(s0)
  lp <- rowSums(log(stats::pnorm((1 - th0) / sd) - stats::pnorm(-th0 / sd)))
  lw <- -n * lp
  expect_equal(res$log_box_prob, lp)
  expect_equal(res$weights, exp(lw - max(lw)) / sum(exp(lw - max(lw))))
  expect_equal(res$n_theta, 3L)
  expect_equal(sum(res$weights), 1)
  expect_equal(res$ess, 1 / sum(res$weights^2))
  expect_equal(res$ess_frac, res$ess / nd)
  expect_equal(res$tilt_nats, n * abs(stats::median(lp)))
  # the hot chain holds Sigma0 = I and theta0 = 0; it must not be read
  expect_equal(length(res$weights), nd)
})


test_that("the general covariance path agrees with the diagonal shortcut", {
  lbp <- impala:::log_box_prob
  mean <- c(0.3, 0.6)
  cov <- diag(c(0.04, 0.09))
  nudged <- cov + matrix(c(0, 1e-9, 1e-9, 0), 2, 2)
  # 1e-9 is above the diagonal tolerance, so this takes the pmvnorm path
  expect_equal(lbp(mean, nudged), lbp(mean, cov), tolerance = 1e-6)

  # and with real correlation it still agrees with a direct pmvnorm call
  skip_if_not_installed("mvtnorm")
  corr <- matrix(c(0.04, 0.03, 0.03, 0.09), 2, 2)
  direct <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(1, 1), mean = mean,
                             sigma = corr)
  expect_equal(lbp(mean, corr), log(as.numeric(direct)), tolerance = 1e-5)

  # p = 1 goes through the univariate branch
  expect_equal(lbp(0.4, matrix(0.09)),
               log(stats::pnorm(0.6 / 0.3) - stats::pnorm(-0.4 / 0.3)))
})


test_that("log_box_prob does not disturb the random number stream", {
  skip_if_not_installed("mvtnorm")
  set.seed(7)
  before <- stats::runif(1)
  seed_after <- .Random.seed
  impala:::log_box_prob(c(0.3, 0.6), matrix(c(0.04, 0.03, 0.03, 0.09), 2, 2))
  expect_identical(.Random.seed, seed_after)
  expect_equal(before, 0.9889093, tolerance = 1e-6)
})


test_that("weights are flat and nothing is corrected when the bounds are slack", {
  setup <- trunc_setup(p = 2, ng = 5)
  set.seed(1)
  nd <- 200
  th0 <- 0.5 + matrix(stats::rnorm(nd * 2, sd = 0.01), nd, 2)
  s0 <- matrix(0.01^2, nd, 2)
  res <- hierTruncationWeights(setup, trunc_out(th0, s0), mcmc_use = 1:nd)
  expect_gt(res$ess_frac, 0.99)
  expect_lt(res$tilt_nats, 1e-6)
  expect_equal(res$theta0_sd_raw, res$theta0_sd_reweighted, tolerance = 1e-3)
  expect_equal(res$theta0_mean_raw, res$theta0_mean_reweighted,
               tolerance = 1e-3)
  expect_equal(res$sigma0_diag_median_raw, res$sigma0_diag_median_reweighted,
               tolerance = 1e-6)
  # raw summaries are the plain sample ones (sd with denominator n)
  expect_equal(res$theta0_mean_raw, colMeans(th0))
  expect_equal(res$theta0_sd_raw,
               sqrt(colMeans(sweep(th0, 2, colMeans(th0))^2)))
})


test_that("reweighting widens Sigma0 when the bounds bind", {
  setup <- trunc_setup(p = 1, ng = 6)
  set.seed(2)
  nd <- 4000
  th0 <- matrix(stats::runif(nd, 0.2, 0.8), nd, 1)
  s0 <- matrix(stats::runif(nd, 0.05, 0.6), nd, 1)^2
  res <- hierTruncationWeights(setup, trunc_out(th0, s0), mcmc_use = 1:nd)
  expect_gt(res$tilt_nats, 0.1)
  # the contraction the correction exists to undo: Sigma0 must go up
  expect_gt(res$sigma0_diag_median_reweighted[1],
            res$sigma0_diag_median_raw[1])
  expect_lt(res$ess_frac, 1)
})


test_that("mcmc_use defaults to the second half of the chain", {
  setup <- trunc_setup(p = 1, ng = 2)
  nd <- 100
  out <- trunc_out(matrix(0.5, nd, 1), matrix(0.01, nd, 1))
  res <- hierTruncationWeights(setup, out)
  expect_equal(length(res$weights), nd / 2)
  expect_equal(length(hierTruncationWeights(setup, out, mcmc_use = 10:20)$weights), 11)
})


test_that("arguments are validated", {
  setup <- trunc_setup(p = 1, ng = 2)
  out <- trunc_out(matrix(0.5, 20, 1), matrix(0.01, 20, 1))
  expect_error(hierTruncationWeights(setup, out, temperature = 3),
               "between 1 and 2")
  expect_error(hierTruncationWeights(setup, out, mcmc_use = 0:5), "indices")
  expect_error(hierTruncationWeights(setup, out, mcmc_use = integer(0)),
               "non-empty")
  expect_error(hierTruncationWeights(setup, list(theta0 = out$theta0)),
               "Sigma0")
  bare <- setup
  bare$theta0_prior_mean <- NULL
  expect_error(hierTruncationWeights(bare, out), "setHierPriors")
  # the hot chain is reachable, and carries Sigma0 = I there
  hot <- hierTruncationWeights(setup, out, mcmc_use = 1:20, temperature = 2)
  expect_equal(hot$log_box_prob,
               rep(log(stats::pnorm(1) - stats::pnorm(0)), 20))
})


test_that("a custom constraint_func warns that the correction is partial", {
  always <- function(x, bounds, constants = NULL) rep(TRUE, length(x[[1]]))
  setup <- trunc_setup(p = 1, ng = 2, constraint_func = always)
  out <- trunc_out(matrix(0.5, 10, 1), matrix(0.04, 10, 1))
  expect_warning(hierTruncationWeights(setup, out, mcmc_use = 1:10),
                 "bounds only")
  # the default constraint function does not warn
  plain <- trunc_setup(p = 1, ng = 2)
  expect_silent(hierTruncationWeights(plain, out, mcmc_use = 1:10))
})


test_that("draws with an unusable box probability are dropped with a warning", {
  setup <- trunc_setup(p = 2, ng = 2)
  nd <- 10
  out <- trunc_out(matrix(0.5, nd, 2), matrix(0.04, nd, 2))
  # a corrupted Sigma0 draw: the box probability cannot be computed
  out$Sigma0[3, 1, , ] <- NA_real_
  expect_warning(res <- hierTruncationWeights(setup, out, mcmc_use = 1:nd),
                 "dropped")
  expect_equal(length(res$weights), nd - 1)
})


test_that("weightedQuantile matches the unweighted sample", {
  x <- c(5, 1, 3, 2, 4)
  w <- rep(0.2, 5)
  expect_equal(weightedQuantile(x, 0.5, w), stats::median(x))
  # equal weights reproduce the type = 5 plotting positions
  expect_equal(weightedQuantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9), w),
               as.numeric(stats::quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9),
                                          type = 5)))
  # duplicating a value by doubling its weight shifts the median toward it
  expect_equal(weightedQuantile(c(1, 2, 3), 0.5, c(1, 1, 1)), 2)
  expect_lt(weightedQuantile(c(1, 2, 3), 0.5, c(5, 1, 1)), 2)
  # clamped to the sample range at the extremes, and NULL weights are equal
  expect_equal(weightedQuantile(x, c(0, 1)), c(1, 5))
  expect_equal(weightedQuantile(x, 0.5), stats::median(x))
  # scale of the weights is irrelevant
  expect_equal(weightedQuantile(x, 0.3, w), weightedQuantile(x, 0.3, 100 * w))
  # zero-weight entries are ignored
  expect_equal(weightedQuantile(c(1, 2, 3, 99), 0.5, c(1, 1, 1, 0)), 2)
  expect_equal(weightedQuantile(c(7, 7), 0.5, c(1, 0)), 7)
  expect_error(weightedQuantile(x, 0.5, 1:2), "same length")
  expect_error(weightedQuantile(x, 1.5, w), "0, 1")
  expect_error(weightedQuantile(x, 0.5, rep(0, 5)), "positive weight")
})


test_that("the diagnostic runs on a calibHier fit", {
  skip_on_cran()
  theta_true <- rbind(c(0.3, 0.6), c(0.5, 0.4), c(0.4, 0.5))
  set.seed(11)
  ng <- 3
  ny_per <- 15
  ny <- ng * ny_per
  A <- matrix(stats::rnorm(ny * 2), ny, 2)
  group <- rep(seq_len(ng), each = ny_per)
  yobs <- rowSums(A * theta_true[group, ]) + stats::rnorm(ny, 0, 0.03)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, yobs, StubModel(A),
                             sd_est = rep(0.03, ng), s2_df = rep(20, ng),
                             s2_ind = group, theta_ind = group)
  setup <- setMCMC(setup, nmcmc = 1500, start_adapt_iter = 200, decor = 50,
                   start_var_theta = 1e-4, start_var_ls2 = 1e-3)
  setup <- setTemperatureLadder(setup, 1.3^(0:1), start_temper = 100)
  setup <- setHierPriors(setup, Sigma0_prior_df = 5,
                         Sigma0_prior_scale = diag(2) * 0.05)
  out <- calibHier(setup)
  res <- hierTruncationWeights(setup, out, mcmc_use = 751:1500)
  expect_equal(length(res$weights), 750)
  expect_equal(res$n_theta, 3L)
  expect_true(all(is.finite(res$log_box_prob)))
  expect_true(all(res$log_box_prob <= 0))
  expect_gt(res$ess_frac, 0)
  expect_lte(res$ess_frac, 1)
  expect_length(res$theta0_mean_reweighted, 2)
  expect_length(res$sigma0_diag_median_reweighted, 2)
  # this fit sits mid-box with a tight Sigma0, so nothing needs correcting
  expect_gt(res$ess_frac, 0.9)
  expect_equal(res$theta0_mean_raw, res$theta0_mean_reweighted,
               tolerance = 0.01)
})


test_that("reweighting recovers the stated hyperprior it was derived from", {
  skip_on_cran()
  # The round trip the correction is for. Draw (theta0, Sigma0) from the stated
  # hyperprior, then resample with probability proportional to P(B)^n -- that
  # is the effective hyperprior calibHier targets, because it applies the
  # truncation without its normalising constant. Reweighting those draws by
  # P(B)^(-n) must give back the summaries of the stated prior. Sigma0 comes
  # from an inverse-Wishart, so the general pmvnorm path is exercised.
  set.seed(3)
  p <- 2
  n <- 6
  N <- 30000
  nd <- 4000
  sd_prior <- 0.2
  th0 <- matrix(0.5 + stats::rnorm(N * p, sd = sd_prior), N, p)
  scale <- diag(p) * 0.3^2
  W <- stats::rWishart(N, p + 4, chol2inv(chol(scale)))
  S0 <- lapply(seq_len(N), function(j) chol2inv(chol(W[, , j])))
  lp <- vapply(seq_len(N),
               function(j) impala:::log_box_prob(th0[j, ], S0[[j]]),
               numeric(1))
  prior_sigma11 <- stats::median(vapply(S0, function(s) s[1, 1], numeric(1)))

  tilt <- exp(n * lp - max(n * lp))
  keep <- sample.int(N, nd, replace = TRUE, prob = tilt / sum(tilt))
  out <- list(theta0 = array(0, c(nd, 1, p)),
              Sigma0 = array(0, c(nd, 1, p, p)))
  out$theta0[, 1, ] <- th0[keep, ]
  for (j in seq_len(nd)) out$Sigma0[j, 1, , ] <- S0[[keep[j]]]

  setup <- trunc_setup(p = 2, ng = n)
  expect_equal(sum(setup$ntheta), n)
  res <- hierTruncationWeights(setup, out, mcmc_use = seq_len(nd))

  # the contraction is real and the correction removes most of it
  expect_true(all(res$theta0_sd_raw < 0.9 * sd_prior))
  expect_equal(res$theta0_sd_reweighted, rep(sd_prior, p), tolerance = 0.06)
  expect_lt(res$sigma0_diag_median_raw[1], 0.92 * prior_sigma11)
  expect_equal(res$sigma0_diag_median_reweighted[1], prior_sigma11,
               tolerance = 0.06)
  # theta0 is centred either way; only its spread is contracted
  expect_equal(res$theta0_mean_reweighted, rep(0.5, p), tolerance = 0.03)
  expect_gt(res$tilt_nats, 0.05)
  expect_gt(res$ess_frac, 0.05)
})
