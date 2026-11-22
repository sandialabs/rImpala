# Ported from python impala's test_friedman_fit.py. The Friedman function uses
# only the first four of six parameters, so the posterior of the other two
# must stay uniform while the four that are used must not.

friedman <- function(grid) {
  # vectorised over rows of theta: returns (nrow x length(grid))
  function(theta) {
    theta <- matrix(theta, ncol = 6)
    10 * sin(pi * outer(theta[, 1], grid)) +
      20 * (theta[, 2] - 0.5)^2 + 10 * theta[, 3] + 5 * theta[, 4]
  }
}

test_that("unused Friedman parameters keep a uniform posterior", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(0)
  n_features <- 6
  grid <- seq(0, 1, length.out = 50)
  f <- friedman(grid)
  theta_true <- stats::runif(n_features, 0.1, 0.9)
  yobs <- as.numeric(f(theta_true)) + stats::rnorm(length(grid), 0, 0.1)

  bounds <- stats::setNames(rep(list(c(0, 1)), n_features), as.character(0:5))
  setup <- CalibSetup(bounds)
  model <- ModelF(f, names(bounds), vectorized = TRUE)
  setup <- addVecExperiments(setup, yobs, model, sd_est = 1, s2_df = 0,
                             s2_ind = rep(1, length(yobs)))
  setup <- setTemperatureLadder(setup, 1.05^(0:19))
  setup <- setMCMC(setup, nmcmc = 15000, decor = 100)
  out <- calibPool(setup)

  post <- out$theta[seq(6001, 15000, by = 5), 1, ]
  # rejected proposals repeat draws, so ks.test warns about ties
  for (j in 1:4) {
    expect_lt(suppressWarnings(stats::ks.test(post[, j], "punif"))$p.value, 0.01)
  }
  for (j in 5:6) {
    expect_gt(suppressWarnings(stats::ks.test(post[, j], "punif"))$p.value, 0.01)
  }
})


test_that("ModelF gives the same predictions vectorised or row by row", {
  grid <- seq(0, 1, length.out = 7)
  f <- friedman(grid)
  parmat <- stats::setNames(lapply(1:6, function(j) stats::runif(3)), as.character(0:5))
  by_row <- evalm(ModelF(function(th) as.numeric(f(th)), names(parmat)), parmat)
  vec <- evalm(ModelF(f, names(parmat), vectorized = TRUE), parmat)
  expect_equal(dim(by_row), c(3, 7))
  expect_equal(by_row, vec)
})
