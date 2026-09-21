# Tests for the discrepancy sampling path -- the nd > 0 branches in
# addVecExperiments() and calibPool() and the StubModel discrep_sample method.
# All tests run without BASS/mvBayes via the linear StubModel from
# helper-stub-model.R, which gains a discrep_sample method there.

# Build a CalibSetup around StubModel with a discrepancy basis D.
stub_setup_discrep <- function(p = 2, ny = 12, ns2 = 1, ntemps = 1, nd = 2,
                               discrep_tau = 1, s2mode = "gibbs", nmcmc = 300,
                               seed = 21) {
  set.seed(seed)
  bounds <- list()
  for (j in 1:p) bounds[[paste0("t_", j)]] <- c(0, 1)

  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- seq(0.3, 0.7, length.out = p)
  # a smooth discrepancy basis (polynomial columns over the index)
  grid <- seq(0, 1, length.out = ny)
  D <- vapply(1:nd, function(k) grid^k, numeric(ny))
  yobs <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, 0.05)

  setup <- CalibSetup(bounds, cf_bounds)
  setup <- addVecExperiments(
    setup, yobs, StubModel(A, s2 = s2mode),
    sd_est = rep(0.05, ns2),
    s2_df  = rep(20, ns2),
    s2_ind = rep(1:ns2, length.out = ny),
    D = D,
    discrep_tau = discrep_tau
  )
  setup <- setMCMC(setup, nmcmc = nmcmc, start_adapt_iter = 100, decor = 50)
  if (ntemps > 1) {
    setup <- setTemperatureLadder(setup, 1.3^(0:(ntemps - 1)),
                                  start_temper = 100)
  }
  list(setup = setup, D = D, nd = nd, nmcmc = nmcmc, ntemps = ntemps, p = p)
}


test_that("addVecExperiments records the discrepancy basis and settings", {
  ny <- 10; nd <- 3
  A <- matrix(stats::rnorm(ny * 2), ny, 2)
  D <- matrix(stats::rnorm(ny * nd), ny, nd)
  mod <- StubModel(A)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, stats::rnorm(ny), mod,
                             sd_est = 0.1, s2_df = 5,
                             s2_ind = rep(1, ny),
                             D = D, discrep_tau = 2.5)

  expect_equal(setup$models[[1]]$nd, nd)
  expect_equal(setup$models[[1]]$D, D)
  expect_equal(setup$models[[1]]$discrep_tau, 2.5)
})


test_that("addVecExperiments stores a measurement error correlation", {
  ny <- 8
  A <- matrix(stats::rnorm(ny * 2), ny, 2)
  mec <- diag(ny); mec[1, 2] <- mec[2, 1] <- 0.3
  mod <- StubModel(A)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, stats::rnorm(ny), mod,
                             sd_est = 0.1, s2_df = 5, s2_ind = rep(1, ny),
                             meas_error_cor = mec)
  expect_equal(setup$models[[1]]$meas_error_cor, mec)
})


test_that("addVecExperiments leaves nd = 0 when no basis is supplied", {
  ny <- 6
  A <- matrix(stats::rnorm(ny * 2), ny, 2)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                             sd_est = 0.1, s2_df = 5, s2_ind = rep(1, ny))
  expect_equal(setup$models[[1]]$nd, 0)
  expect_null(setup$models[[1]]$D)
})


test_that("calibPool records discrepancy coefficients of the right shape", {
  s <- stub_setup_discrep(ntemps = 1, nd = 2)
  set.seed(101)
  out <- suppressMessages(calibPool(s$setup))

  # discrep_vars is (nmcmc x ntemps x nd), one experiment
  expect_length(out$discrep_vars, 1)
  expect_equal(dim(out$discrep_vars[[1]]), c(s$nmcmc, s$ntemps, s$nd))

  # the first row is the initialised zero; later draws must not all be zero
  post <- out$discrep_vars[[1]][-1, 1, , drop = FALSE]
  expect_true(any(post != 0))
  expect_false(any(is.na(post)))
  expect_true(all(is.finite(post)))
})


test_that("calibPool runs the discrepancy branch under tempering", {
  s <- stub_setup_discrep(ntemps = 3, nd = 2, nmcmc = 300)
  set.seed(202)
  out <- suppressMessages(calibPool(s$setup))

  expect_equal(dim(out$discrep_vars[[1]]), c(s$nmcmc, 3, s$nd))
  # every temperature should have moved off the zero start at some point
  for (t in 1:3) {
    expect_true(any(out$discrep_vars[[1]][-1, t, ] != 0),
                info = paste("temperature", t))
  }
  # a swap having occurred is not required, but the run must stay finite
  expect_true(all(is.finite(out$theta)))
  expect_true(all(is.finite(out$llik)))
})


test_that("calibPool discrepancy path works with an M-H s2 update", {
  s <- stub_setup_discrep(ntemps = 1, nd = 2, s2mode = "mh", nmcmc = 300)
  set.seed(303)
  out <- suppressMessages(calibPool(s$setup))
  expect_equal(dim(out$discrep_vars[[1]]), c(s$nmcmc, 1, s$nd))
  expect_true(all(is.finite(out$s2[[1]])))
})


test_that("discrep_sample returns a finite length-nd coefficient vector", {
  ny <- 10; nd <- 3
  A <- matrix(stats::rnorm(ny * 2), ny, 2)
  D <- vapply(1:nd, function(k) seq(0, 1, length.out = ny)^k, numeric(ny))
  mod <- StubModel(A)
  mod$D <- D; mod$nd <- nd; mod$discrep_tau <- 1

  cov <- lik_cov_inv.StubModel(mod, rep(0.05^2, ny))
  yobs <- stats::rnorm(ny)
  pred <- stats::rnorm(ny)

  set.seed(7)
  b <- discrep_sample.StubModel(mod, yobs, pred, cov, itl = 1)
  expect_length(b, nd)
  expect_true(all(is.finite(b)))
})
