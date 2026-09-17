# fixTheta() pins one or more calibration parameters so calibPool samples the
# conditional posterior of the rest -- the mechanism behind conditional and
# cut-Bayes inference. Two things must hold: the fixed columns never move, and
# the free ones target the same conditional posterior a smaller model would.

test_that("fixTheta validates its arguments", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(-2, 2)), cf_bounds)

  expect_error(fixTheta(setup), "character vector")
  expect_error(fixTheta(setup, 1, 0.5), "character vector")
  expect_error(fixTheta(setup, "nope", 0.5), "not in set of input names")
  expect_error(fixTheta(setup, c("t_1", "t_1"), c(0.2, 0.3)), "duplicate")
  expect_error(fixTheta(setup, "t_1", c(0.2, 0.3)), "same length as pname")
  expect_error(fixTheta(setup, "t_1", "0.5"), "numeric vector")
  expect_error(fixTheta(setup, "t_1", NA_real_), "must be finite")
  # bounds are checked on the native scale, per parameter
  expect_error(fixTheta(setup, "t_1", 1.5), "outside of bounds")
  expect_error(fixTheta(setup, "t_2", -3), "outside of bounds")
  expect_silent(fixTheta(setup, "t_2", -1.5))
})


test_that("fixTheta records values in bounds order and replaces on re-fix", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1), t_3 = c(0, 1)),
                      cf_bounds)

  setup <- fixTheta(setup, "t_3", 0.9)
  expect_equal(setup$theta_fixed, c(t_3 = 0.9))

  # added out of order, stored in bounds order
  setup <- fixTheta(setup, "t_1", 0.1)
  expect_named(setup$theta_fixed, c("t_1", "t_3"))

  # re-fixing replaces rather than duplicating
  setup <- fixTheta(setup, "t_1", 0.2)
  expect_equal(setup$theta_fixed, c(t_1 = 0.2, t_3 = 0.9))
})


test_that("fixTheta warns about a prior on a fixed parameter", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addThetaPrior(setup, "normal", list(mean = 0.5, sd = 0.2), "t_1")

  expect_warning(fixTheta(setup, "t_1", 0.5), "Prior set for fixed parameter")
  # a prior on the other parameter is fine
  expect_silent(fixTheta(setup, "t_2", 0.5))

  # joint priors are checked through their `names` field too
  joint <- addJointThetaPrior(setup, c("t_1", "t_2"), function(p) p$t_1)
  expect_warning(fixTheta(joint, "t_2", 0.5), "Prior set for fixed parameter")
})


test_that("fixTheta warns when every parameter is fixed", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  expect_warning(fixTheta(setup, c("t_1", "t_2"), c(0.3, 0.4)),
                 "All 2 calibration parameters are fixed")
})


test_that("unfixTheta releases parameters", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- fixTheta(setup, "t_1", 0.3)
  # fixing the second of two parameters warns that nothing is left to sample
  setup <- suppressWarnings(fixTheta(setup, "t_2", 0.4))

  one <- unfixTheta(setup, "t_1")
  expect_equal(one$theta_fixed, c(t_2 = 0.4))

  # releasing the last one clears the field entirely
  expect_null(unfixTheta(one, "t_2")$theta_fixed)
  # no argument releases everything
  expect_null(unfixTheta(setup)$theta_fixed)

  expect_error(unfixTheta(setup, "nope"), "not in set of input names")
  expect_warning(unfixTheta(unfixTheta(setup), "t_1"), "not currently fixed")
})


test_that("fixed_theta_split maps native values onto the unit scale", {
  fixed_theta_split <- getFromNamespace("fixed_theta_split", "impala")
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(-2, 2), t_3 = c(10, 20)),
                      cf_bounds)

  # nothing fixed: every parameter free, no unit values
  s0 <- fixed_theta_split(setup)
  expect_equal(s0$free_idx, 1:3)
  expect_length(s0$fixed_idx, 0)

  setup <- fixTheta(setup, c("t_2", "t_3"), c(0, 15))
  s <- fixed_theta_split(setup)
  expect_equal(s$fixed_idx, c(2L, 3L))
  expect_equal(s$free_idx, 1L)
  # (0 - -2)/4 = 0.5 and (15 - 10)/10 = 0.5
  expect_equal(as.numeric(s$unit), c(0.5, 0.5))
})


test_that("fixed_theta_split treats a legacy setup as having nothing fixed", {
  # setups pickled before theta_fixed existed have no such field at all
  fixed_theta_split <- getFromNamespace("fixed_theta_split", "impala")
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup$theta_fixed <- NULL
  expect_equal(fixed_theta_split(setup)$free_idx, 1:2)
})


test_that("update_m/gen_cand adapt only the requested columns", {
  set.seed(51)
  ntemps <- 3; p <- 4; free <- c(2L, 4L); nmcmc <- 40
  x <- array(stats::runif(nmcmc * ntemps * p), dim = c(nmcmc, ntemps, p))
  # fixed columns are constant in a real run; make them so here
  x[, , c(1, 3)] <- 0.5

  obj <- getFromNamespace("AMcov_pool", "impala")(ntemps, length(free),
                                                 1e-4, 10, 0)
  obj <- update_m(obj, x, 10, cols = free)
  expect_equal(dim(obj$cov), c(ntemps, length(free), length(free)))
  # the running mean must come from the free columns, not the first two
  for (t in 1:ntemps) {
    expect_equal(obj$mu[t, ], colMeans(x[1:10, t, free]), tolerance = 1e-10)
  }

  cand <- gen_cand(obj, x, 12, cols = free)
  expect_equal(dim(cand), c(ntemps, length(free)))

  # candidates are centred on the previous free values, not on the constant
  # columns; use an unadapted, very tight proposal so the centre is unambiguous
  tight <- getFromNamespace("AMcov_pool", "impala")(ntemps, length(free),
                                                   1e-12, 1000, 0)
  cand <- gen_cand(tight, x, 12, cols = free)
  expect_equal(cand, matrix(x[11, , free], ntemps, length(free)),
               tolerance = 1e-4)
  # and would differ from the fixed columns, which sit at 0.5
  expect_gt(max(abs(cand - 0.5)), 1e-3)
})


test_that("calibPool holds fixed parameters constant and samples the rest", {
  set.seed(61)
  p <- 3; ny <- 40; sig <- 0.05
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.3, 0.5, 0.7)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)

  bounds <- list(t_1 = c(0, 1), t_2 = c(0, 1), t_3 = c(0, 1))
  setup <- CalibSetup(bounds, cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                             s2_df = 20000, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 3000, start_adapt_iter = 300, decor = 100)
  # fix the middle parameter at its true value
  setup <- fixTheta(setup, "t_2", theta_true[2])

  out <- suppressWarnings(calibPool(setup))

  expect_equal(dim(out$theta), c(3000, 1, p))
  expect_equal(out$theta_fixed, c(t_2 = theta_true[2]))

  # the fixed column never moves, in either scale, at any iteration
  expect_true(all(out$theta[, 1, 2] == theta_true[2]))
  expect_true(all(out$theta_native$t_2 == theta_true[2]))

  # the free columns still mix and recover the truth
  keep <- 1501:3000
  for (k in c(1, 3)) {
    expect_gt(stats::sd(out$theta[keep, 1, k]), 0)
    expect_lt(abs(mean(out$theta[keep, 1, k]) - theta_true[k]), 0.1)
  }

  # the decorrelation step must not touch the fixed parameter either
  expect_equal(out$count_decor[2, ], 0)
  expect_gt(sum(out$count_decor[c(1, 3), ]), 0)
})


test_that("fixing a parameter matches calibrating a reduced model", {
  # Fixing t_2 at v should give exactly the posterior of the p=2 problem with
  # A[, 2] * v folded into the data. This is the real correctness check: the
  # conditional posterior, not just a frozen column.
  set.seed(67)
  p <- 3; ny <- 60; sig <- 0.05
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.35, 0.60, 0.45)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)
  v <- 0.6   # value t_2 is held at

  # analytic conditional posterior for the two free parameters, flat prior
  Af <- A[, c(1, 3)]
  yf <- y - A[, 2] * v
  AtA_inv <- solve(crossprod(Af))
  mu_an <- as.numeric(AtA_inv %*% crossprod(Af, yf))
  sd_an <- sqrt(diag(sig^2 * AtA_inv))

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1), t_3 = c(0, 1)),
                      cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                             s2_df = 20000, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 6000, start_adapt_iter = 500, decor = 100)
  setup <- fixTheta(setup, "t_2", v)

  out <- suppressWarnings(calibPool(setup))
  keep <- 3001:6000
  th <- out$theta[keep, 1, c(1, 3)]

  # same tolerances as test-posterior-correctness.R: mean within half an
  # analytic sd, spread within 30%
  expect_lt(max(abs(colMeans(th) - mu_an) / sd_an), 0.5)
  expect_lt(max(abs(apply(th, 2, stats::sd) / sd_an - 1)), 0.3)
})


test_that("fixing every parameter still samples s2", {
  set.seed(71)
  p <- 2; ny <- 200; sig <- 0.08
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.4, 0.6)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = 0.2,
                             s2_df = 2, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 2000, start_adapt_iter = 300, decor = 100)
  setup <- suppressWarnings(fixTheta(setup, c("t_1", "t_2"), theta_true))

  out <- suppressWarnings(calibPool(setup))

  # theta is pinned everywhere
  expect_true(all(out$theta[, 1, 1] == theta_true[1]))
  expect_true(all(out$theta[, 1, 2] == theta_true[2]))
  expect_equal(sum(out$count), 0)
  expect_equal(sum(out$count_decor), 0)

  # ... but s2 still travels from the deliberately wrong 0.2 down to the truth
  keep <- 1001:2000
  expect_lt(abs(sqrt(mean(out$s2[[1]][keep, 1, ])) - sig) / sig, 0.4)
  expect_true(all(is.finite(out$llik[-1])))
})


test_that("fixed parameters survive tempering swaps", {
  set.seed(73)
  p <- 3; ny <- 30; sig <- 0.05
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.3, 0.5, 0.7)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1), t_3 = c(0, 1)),
                      cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                             s2_df = 20000, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 1500, start_adapt_iter = 200, decor = 50)
  setup <- setTemperatureLadder(setup, 1.3^(0:3), start_temper = 300)
  setup <- fixTheta(setup, c("t_1", "t_3"), c(0.3, 0.7))

  out <- suppressWarnings(calibPool(setup))

  # a swap moves whole theta rows between temperatures; because every chain
  # shares the same fixed values the pinned columns stay pinned at all temps
  expect_true(all(out$theta[, , 1] == 0.3))
  expect_true(all(out$theta[, , 3] == 0.7))
  expect_gt(stats::sd(out$theta[, 1, 2]), 0)
  expect_lt(abs(mean(out$theta[751:1500, 1, 2]) - theta_true[2]), 0.15)
})


test_that("a single free parameter out of many works", {
  # p_free == 1 exercises the degenerate proposal shape, as p == 1 does normally
  set.seed(79)
  p <- 3; ny <- 30; sig <- 0.05
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.3, 0.5, 0.7)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1), t_3 = c(0, 1)),
                      cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                             s2_df = 20000, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 1500, start_adapt_iter = 200, decor = 50)
  setup <- fixTheta(setup, c("t_1", "t_2"), theta_true[1:2])

  out <- suppressWarnings(calibPool(setup))
  expect_equal(dim(out$cov_theta_cand$S), c(1, 1, 1))
  expect_lt(abs(mean(out$theta[751:1500, 1, 3]) - theta_true[3]), 0.1)
})


test_that("unfixTheta restores full sampling in calibPool", {
  set.seed(83)
  fx <- stub_setup(3, 20, 1, 1, nmcmc = 800)
  setup <- fixTheta(fx$setup, "t_2", 0.5)
  setup <- unfixTheta(setup, "t_2")

  out <- suppressWarnings(calibPool(setup))
  expect_null(out$theta_fixed)
  # every column moves again
  for (k in 1:3) expect_gt(stats::sd(out$theta[401:800, 1, k]), 0)
})
