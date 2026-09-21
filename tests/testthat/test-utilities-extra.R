# Additional unit tests for internal helpers not covered by test-utilities.R:
# the full set of eval_theta_priors distributions, the Sherman-Woodbury-Morrison
# identity (swm), the Gaussian draw helper (chol_sample), the bounded
# inverse-gamma sampler's retry/clamp logic (rig_bounded), ndims, and the
# single-row-bounds branches of unnormalize.

test_that("ndims returns the number of array dimensions", {
  ndims <- getFromNamespace("ndims", "impala")
  expect_equal(ndims(1:5), 0L)                         # a bare vector has no dim
  expect_equal(ndims(matrix(0, 2, 3)), 2L)
  expect_equal(ndims(array(0, dim = c(2, 3, 4))), 3L)
})


test_that("eval_theta_priors evaluates every supported distribution", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")

  x <- c(0.2, 0.5, 0.8)
  theta <- list(t_1 = x)

  cases <- list(
    list(dist = "normal",
         params = list(mean = 0.5, sd = 0.3),
         ref = stats::dnorm(x, 0.5, 0.3, log = TRUE)),
    list(dist = "lognormal",
         params = list(meanlog = -1, sdlog = 0.5),
         ref = stats::dlnorm(x, -1, 0.5, log = TRUE)),
    list(dist = "beta",
         params = list(shape1 = 2, shape2 = 5),
         ref = stats::dbeta(x, 2, 5, log = TRUE)),
    list(dist = "uniform",
         params = list(min = 0, max = 1),
         ref = stats::dunif(x, 0, 1, log = TRUE)),
    list(dist = "gamma",
         params = list(shape = 2, rate = 3),
         ref = stats::dgamma(x, shape = 2, rate = 3, log = TRUE)),
    list(dist = "cauchy",
         params = list(location = 0.5, scale = 0.2),
         ref = stats::dcauchy(x, 0.5, 0.2, log = TRUE))
  )

  for (case in cases) {
    priors <- list(list(name = "t_1", dist = case$dist, params = case$params))
    got <- eval_theta_priors(theta, priors)
    expect_equal(got, case$ref, tolerance = 1e-12,
                 info = paste("distribution:", case$dist))
  }
})


test_that("eval_theta_priors sums independent priors and returns zeros with none", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")

  x1 <- c(0.3, 0.6); x2 <- c(0.4, 0.9)
  theta <- list(t_1 = x1, t_2 = x2)

  # no priors -> a length-matching vector of zeros
  expect_equal(eval_theta_priors(theta, list()), c(0, 0))
  expect_equal(eval_theta_priors(theta, NULL), c(0, 0))

  priors <- list(
    list(name = "t_1", dist = "normal", params = list(mean = 0.5, sd = 0.3)),
    list(name = "t_2", dist = "beta",   params = list(shape1 = 2, shape2 = 2))
  )
  got <- eval_theta_priors(theta, priors)
  want <- stats::dnorm(x1, 0.5, 0.3, log = TRUE) +
          stats::dbeta(x2, 2, 2, log = TRUE)
  expect_equal(got, want, tolerance = 1e-12)
})


test_that("eval_theta_priors honours the tnames renaming argument", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")

  x <- c(0.25, 0.75)
  # theta is supplied without names; tnames assigns them so the prior can find
  # its parameter by name
  theta <- list(x)
  priors <- list(list(name = "alpha", dist = "normal",
                      params = list(mean = 0.5, sd = 0.4)))

  got <- eval_theta_priors(theta, priors, tnames = "alpha")
  expect_equal(got, stats::dnorm(x, 0.5, 0.4, log = TRUE), tolerance = 1e-12)
})


test_that("eval_theta_priors errors on an unsupported distribution", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")
  theta <- list(t_1 = c(0.5))
  priors <- list(list(name = "t_1", dist = "weibull", params = list()))
  expect_error(eval_theta_priors(theta, priors), "Unsupported dist")
})


test_that("eval_theta_priors dispatches to a joint prior's density function", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")

  theta <- list(t_1 = c(0.3, 0.6), t_2 = c(0.4, 0.1))
  # a joint prior is recognised by its `names` field and calls log_density_fn
  joint <- list(list(
    names = c("t_1", "t_2"),
    log_density_fn = function(p) -0.5 * (p$t_1^2 + p$t_2^2)
  ))
  got <- eval_theta_priors(theta, joint)
  want <- -0.5 * (theta$t_1^2 + theta$t_2^2)
  expect_equal(got, want, tolerance = 1e-12)
})


test_that("swm matches the explicit inverse and log-determinant", {
  swm <- getFromNamespace("swm", "impala")
  chol_solve <- getFromNamespace("chol_solve", "impala")

  set.seed(11)
  n <- 5; k <- 2
  # Build A (n x n) SPD and C (k x k) SPD, with U = t(V) so that
  # M = A + U C V is symmetric positive definite.
  A <- crossprod(matrix(stats::rnorm(n * n), n)) + diag(n)
  C <- crossprod(matrix(stats::rnorm(k * k), k)) + diag(k)
  U <- matrix(stats::rnorm(n * k), n, k)
  V <- t(U)

  Ainv <- solve(A); Cinv <- solve(C)
  Aldet <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
  Cldet <- as.numeric(determinant(C, logarithm = TRUE)$modulus)

  got <- swm(Ainv, U, Cinv, V, Aldet, Cldet)

  M <- A + U %*% C %*% V
  expect_equal(got$inv, solve(M), tolerance = 1e-8)
  expect_equal(got$ldet,
               as.numeric(determinant(M, logarithm = TRUE)$modulus),
               tolerance = 1e-8)
  expect_equal(got$inv %*% M, diag(n), tolerance = 1e-8)
})


test_that("chol_sample recovers the target mean and covariance", {
  chol_sample <- getFromNamespace("chol_sample", "impala")

  set.seed(3)
  mean <- c(1, -2, 0.5)
  cov <- crossprod(matrix(stats::rnorm(9), 3)) + diag(3)

  draws <- replicate(20000, as.numeric(chol_sample(mean, cov)))
  expect_equal(dim(draws), c(3, 20000))
  expect_equal(rowMeans(draws), mean, tolerance = 0.05)
  expect_equal(stats::cov(t(draws)), cov, tolerance = 0.1)
})


test_that("rig_bounded returns draws inside the variance bounds", {
  rig_bounded <- getFromNamespace("rig_bounded", "impala")

  set.seed(9)
  n <- 200
  shape <- rep(3, n); scale <- rep(1, n)
  sd_lower <- rep(0.2, n); sd_upper <- rep(2.0, n)

  s2 <- rig_bounded(shape, scale, sd_lower, sd_upper)
  expect_length(s2, n)
  expect_true(all(s2 >= sd_lower^2 - 1e-12))
  expect_true(all(s2 <= sd_upper^2 + 1e-12))
})


test_that("rig_bounded clamps to the bounds when maxit is exhausted", {
  rig_bounded <- getFromNamespace("rig_bounded", "impala")

  # An impossibly tight window forces the retry loop to exhaust maxit and fall
  # through to the clamp: every returned value must be pinned to a bound.
  set.seed(1)
  n <- 50
  lo <- 0.999; hi <- 1.001
  s2 <- rig_bounded(shape = rep(2, n), scale = rep(1, n),
                    sd_lower = rep(lo, n), sd_upper = rep(hi, n),
                    maxit = 1)
  expect_true(all(s2 >= lo^2 - 1e-12 & s2 <= hi^2 + 1e-12))
  # with maxit = 0 the loop never runs, so anything out of range is clamped
  # rather than resampled
  s2b <- rig_bounded(shape = rep(2, n), scale = rep(1, n),
                     sd_lower = rep(lo, n), sd_upper = rep(hi, n),
                     maxit = 0)
  expect_true(all(s2b >= lo^2 - 1e-12 & s2b <= hi^2 + 1e-12))
})


test_that("rig_bounded is unconstrained with default 0/Inf bounds", {
  rig_bounded <- getFromNamespace("rig_bounded", "impala")

  set.seed(2)
  n <- 10
  s2 <- rig_bounded(shape = rep(2, n), scale = rep(1, n),
                    sd_lower = rep(0, n), sd_upper = rep(Inf, n))
  # a plain inverse-gamma draw with the same seed should be reproduced exactly:
  # no rejection or clamping happens when the bounds are vacuous
  set.seed(2)
  ref <- 1 / stats::rgamma(n, shape = rep(2, n), scale = rep(1, n))
  expect_equal(s2, ref, tolerance = 1e-12)
})


test_that("unnormalize handles a single row of z with multi-row bounds", {
  normalize <- getFromNamespace("normalize", "impala")
  unnormalize <- getFromNamespace("unnormalize", "impala")

  bounds <- rbind(c(-2, 5), c(10, 20), c(0, 1))
  z <- matrix(c(0, 0.5, 1), nrow = 1)   # one point in the unit cube

  x <- unnormalize(z, bounds)
  expect_equal(dim(x), c(1, 3))
  expect_equal(as.numeric(x), c(-2, 15, 1), tolerance = 1e-12)
  expect_equal(normalize(x, bounds), z, tolerance = 1e-12)
})
