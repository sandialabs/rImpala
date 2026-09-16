# Shape tests confirm calibPool runs; these confirm it targets the *right*
# posterior. For the linear-Gaussian StubModel with a flat prior the posterior is
# available in closed form, so the MCMC output can be checked against it directly
# rather than against a previously recorded value.
#
# Iteration counts are set so the Monte Carlo error sits well inside the
# tolerances below (measured max|z| ~ 0.07 against a 0.5 limit), while keeping
# the suite quick enough to run on every CI platform. Still skipped on CRAN.

test_that("theta posterior matches the analytic conjugate posterior", {
  skip_on_cran()

  set.seed(7)
  p <- 2; ny <- 60; sig <- 0.05
  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.45, 0.55)
  y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sig)

  # flat-prior Gaussian posterior for known sigma
  AtA_inv <- solve(crossprod(A))
  mu_an <- as.numeric(AtA_inv %*% crossprod(A, y))
  cov_an <- sig^2 * AtA_inv

  bounds <- list(t_1 = c(0, 1), t_2 = c(0, 1))
  setup <- CalibSetup(bounds, cf_bounds)
  # a very tight s2 prior at the true variance approximates known-sigma
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                             s2_df = 20000, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 6000, start_adapt_iter = 500, decor = 100)

  out <- suppressWarnings(calibPool(setup))
  keep <- 3001:6000
  th <- matrix(out$theta[keep, 1, ], length(keep), p)

  mu_mc <- colMeans(th)
  cov_mc <- stats::cov(th)
  sd_an <- sqrt(diag(cov_an))

  # posterior mean within a fraction of an analytic sd
  expect_lt(max(abs(mu_mc - mu_an) / sd_an), 0.5)
  # posterior spread within 30% of analytic
  expect_lt(max(abs(sqrt(diag(cov_mc)) / sd_an - 1)), 0.3)
  # and the parameter correlation is reproduced
  corr_an <- cov_an[1, 2] / prod(sd_an)
  corr_mc <- cov_mc[1, 2] / prod(sqrt(diag(cov_mc)))
  expect_lt(abs(corr_mc - corr_an), 0.1)
})


test_that("s2 is recovered under a weak prior", {
  skip_on_cran()

  set.seed(19)
  p <- 2; ny <- 60; sig <- 0.08
  A <- matrix(stats::rnorm(ny * p), ny, p)
  y <- as.numeric(A %*% c(0.45, 0.55)) + stats::rnorm(ny, 0, sig)

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  # start deliberately away from the truth so this tests the update, not the init
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = 0.2,
                             s2_df = 2, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = 6000, start_adapt_iter = 500, decor = 100)

  out <- suppressWarnings(calibPool(setup))
  keep <- 3001:6000
  sd_post <- sqrt(mean(out$s2[[1]][keep, 1, ]))
  expect_lt(abs(sd_post - sig) / sig, 0.25)
})


test_that("each error group recovers its own variance", {
  skip_on_cran()

  set.seed(23)
  p <- 2; ny <- 60
  A <- matrix(stats::rnorm(ny * p), ny, p)
  sd_grp <- c(0.03, 0.15)
  s2_ind <- rep(1:2, each = ny / 2)
  y <- as.numeric(A %*% c(0.45, 0.55)) + stats::rnorm(ny, 0, sd_grp[s2_ind])

  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A), sd_est = c(0.1, 0.1),
                             s2_df = c(2, 2), s2_ind = s2_ind)
  setup <- setMCMC(setup, nmcmc = 6000, start_adapt_iter = 500, decor = 100)

  out <- suppressWarnings(calibPool(setup))
  keep <- 3001:6000
  sd_post <- sqrt(colMeans(out$s2[[1]][keep, 1, ]))

  # both groups distinguished, despite starting from a common 0.1
  expect_lt(max(abs(sd_post - sd_grp) / sd_grp), 0.35)
  expect_lt(sd_post[1], sd_post[2])
})


test_that("tempered and untempered runs agree on the cold-chain posterior", {
  skip_on_cran()

  run <- function(ntemps) {
    set.seed(31)
    p <- 2; ny <- 40; sig <- 0.05
    A <- matrix(stats::rnorm(ny * p), ny, p)
    y <- as.numeric(A %*% c(0.45, 0.55)) + stats::rnorm(ny, 0, sig)
    setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
    setup <- addVecExperiments(setup, y, StubModel(A), sd_est = sig,
                               s2_df = 20000, s2_ind = rep(1, ny))
    setup <- setMCMC(setup, nmcmc = 8000, start_adapt_iter = 500, decor = 100)
    if (ntemps > 1) {
      setup <- setTemperatureLadder(setup, 1.2^(0:(ntemps - 1)),
                                    start_temper = 600)
    }
    out <- suppressWarnings(calibPool(setup))
    keep <- 4001:8000
    colMeans(matrix(out$theta[keep, 1, ], length(keep), p))
  }

  # tempering must not shift the target distribution of the cold chain
  expect_lt(max(abs(run(1) - run(4))), 0.02)
})
