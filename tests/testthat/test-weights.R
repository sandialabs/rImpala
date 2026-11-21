# Observation weights (`wt` in addVecExperiments) raise each observation's
# likelihood contribution to the power of its weight.

weighted_setup <- function(wt, s2mode = "fix", nmcmc = 3000, seed = 7,
                           A = NULL, yobs = NULL, supports_wt = TRUE) {
  set.seed(seed)
  p <- 2
  ny <- 30
  if (is.null(A)) A <- matrix(stats::rnorm(ny * p), ny, p)
  if (is.null(yobs)) yobs <- as.numeric(A %*% c(0.4, 0.6)) + stats::rnorm(nrow(A), 0, 0.1)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, yobs, StubModel(A, s2 = s2mode, supports_wt = supports_wt),
                             sd_est = 0.1, s2_df = 20, s2_ind = rep(1, length(yobs)),
                             wt = wt)
  setup <- setMCMC(setup, nmcmc = nmcmc, start_adapt_iter = 200, decor = 100,
                   start_var_theta = 1e-4)
  list(setup = setup, A = A, yobs = yobs)
}


test_that("addVecExperiments validates and stores weights", {
  A <- matrix(1, 4, 2)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  s1 <- addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4))
  expect_equal(s1$wt[[1]], rep(1, 4))
  expect_error(addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4), wt = 1:3),
               "one weight per element")
  expect_error(addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4), wt = c(1, 1, -1, 1)),
               "non-negative")
})


test_that("non-unit weights need a model that applies them", {
  ws <- weighted_setup(rep(2, 30), supports_wt = FALSE, nmcmc = 10)
  expect_error(calibPool(ws$setup), "does not apply them")
})


test_that("unit weights reproduce the unweighted sampler exactly", {
  a <- weighted_setup(NULL, s2mode = "gibbs", nmcmc = 300)
  b <- weighted_setup(rep(1, 30), s2mode = "gibbs", nmcmc = 300, supports_wt = FALSE)
  set.seed(99)
  out_a <- calibPool(a$setup)
  set.seed(99)
  out_b <- calibPool(b$setup)
  expect_identical(out_a$theta, out_b$theta)
  expect_identical(out_a$s2, out_b$s2)
})


test_that("weighted posterior matches the closed form (fixed s2)", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  wt <- rep(c(3, 0.5), each = 15)
  ws <- weighted_setup(wt, s2mode = "fix", nmcmc = 6000)
  out <- calibPool(ws$setup)
  keep <- 2001:6000
  draws <- out$theta[keep, 1, ]

  # flat prior, Gaussian likelihood raised to wt: precision A' W A / s2
  s2 <- 0.1^2
  prec <- crossprod(ws$A, wt * ws$A) / s2
  post_cov <- solve(prec)
  post_mean <- post_cov %*% crossprod(ws$A, wt * ws$yobs) / s2

  expect_equal(colMeans(draws), as.numeric(post_mean), tolerance = 0.02)
  expect_equal(apply(draws, 2, stats::sd), sqrt(diag(post_cov)), tolerance = 0.25)
})


test_that("weight 2 matches duplicating the data, including the Gibbs s2 draw", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(3)
  A <- matrix(stats::rnorm(20 * 2), 20, 2)
  y <- as.numeric(A %*% c(0.4, 0.6)) + stats::rnorm(20, 0, 0.1)

  w <- weighted_setup(rep(2, 20), s2mode = "gibbs", nmcmc = 5000, A = A, yobs = y)
  d <- weighted_setup(NULL, s2mode = "gibbs", nmcmc = 5000,
                      A = rbind(A, A), yobs = c(y, y))
  set.seed(11)
  out_w <- calibPool(w$setup)
  set.seed(12)
  out_d <- calibPool(d$setup)
  keep <- 1001:5000

  expect_equal(colMeans(out_w$theta[keep, 1, ]), colMeans(out_d$theta[keep, 1, ]), tolerance = 0.02)
  expect_equal(mean(out_w$s2[[1]][keep, 1, 1]), mean(out_d$s2[[1]][keep, 1, 1]), tolerance = 0.1)
})


test_that("ModelBassPca_func applies weights and observation subsets", {
  set.seed(21)
  ny <- 8
  npc <- 2
  obj <- structure(list(
    basis = matrix(stats::rnorm(ny * npc), ny, npc),
    trunc_error_var = seq(1e-3, 2e-3, length.out = ny),
    emu_vars = c(0.02, 0.01),
    npc = npc
  ), class = "ModelBassPca_func")
  s2 <- seq(0.01, 0.03, length.out = ny)
  wt <- c(1, 2, 0.5, 1, 3, 1, 1, 0.25)

  full <- lik_cov_inv(obj, s2)
  mat <- diag(obj$trunc_error_var + s2) + obj$basis %*% diag(obj$emu_vars) %*% t(obj$basis)
  expect_equal(full$inv, solve(mat), tolerance = 1e-8)
  expect_equal(full$ldet, as.numeric(determinant(mat)$modulus), tolerance = 1e-8)

  # weights scale only the diagonal part of the log determinant
  wcov <- lik_cov_inv(obj, s2, wt = wt)
  expect_equal(wcov$inv, full$inv)
  expect_equal(wcov$ldet - full$ldet, sum((wt - 1) * log(obj$trunc_error_var + s2)), tolerance = 1e-10)

  r <- stats::rnorm(ny)
  ll <- llik(obj, r, rep(0, ny), wcov, wt = wt)
  v <- sqrt(wt) * r
  expect_equal(as.numeric(ll), -0.5 * (wcov$ldet + sum(v * (wcov$inv %*% v))))

  # a subset of observations is the covariance of that subset
  inds <- c(2, 5, 7)
  sub <- lik_cov_inv(obj, s2[inds], inds = inds)
  expect_equal(sub$inv, solve(mat[inds, inds]), tolerance = 1e-8)
})


test_that("mvBayes covariances honour an observation subset", {
  set.seed(4)
  ny <- 6
  obj <- structure(list(
    meas_error_corr = diag(ny),
    trunc_error_var = diag(ny) * 1e-4,
    discrep_cov = diag(ny) * 1e-12,
    basis = matrix(stats::rnorm(ny * 2), ny, 2),
    emu_vars = c(1e-3, 2e-3),
    npc = 2
  ), class = "ModelmvBayes")
  s2 <- seq(0.01, 0.02, length.out = ny)
  full <- lik_cov_inv(obj, s2)
  inds <- c(1, 4, 6)
  sub <- lik_cov_inv(obj, s2[inds], inds = inds)
  mat <- solve(full$inv)
  expect_equal(sub$inv, solve(mat[inds, inds]), tolerance = 1e-6)
})
