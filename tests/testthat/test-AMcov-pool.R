# The adaptive Metropolis covariance carries (ntemps x p x p) state. Single-index
# slices of the sample array drop dimensions in R, which used to break every
# ntemps == 1 or p == 1 run, so the degenerate shapes are covered explicitly.

make_cov <- function(ntemps, p, start_var = 1e-4, start_adapt = 10, tau = 0) {
  getFromNamespace("AMcov_pool", "impala")(ntemps, p, start_var, start_adapt, tau)
}

test_that("AMcov_pool initialises with the right shapes", {
  obj <- make_cov(3, 2)
  expect_s3_class(obj, "AMcov_pool")
  expect_equal(dim(obj$S), c(3, 2, 2))
  expect_equal(dim(obj$cov), c(3, 2, 2))
  expect_equal(dim(obj$mu), c(3, 2))
  expect_length(obj$count_100, 3)
  for (t in 1:3) expect_equal(obj$S[t, , ], diag(2) * 1e-4)
})


for (cfg in list(c(1, 1), c(1, 3), c(4, 1), c(4, 3))) {
  ntemps <- cfg[1]; p <- cfg[2]
  test_that(sprintf("update_m/gen_cand keep shapes for ntemps=%d p=%d", ntemps, p), {
    set.seed(100 + 10 * ntemps + p)
    nmcmc <- 40
    x <- array(stats::runif(nmcmc * ntemps * p), dim = c(nmcmc, ntemps, p))
    obj <- make_cov(ntemps, p, start_adapt = 10)

    # m == start_adapt_iter takes the colMeans/cov_3d_pcm branch
    obj <- update_m(obj, x, 10)
    expect_equal(dim(obj$mu), c(ntemps, p))
    expect_equal(dim(obj$cov), c(ntemps, p, p))
    expect_equal(dim(obj$S), c(ntemps, p, p))
    expect_true(all(is.finite(obj$S)))

    # m > start_adapt_iter takes the running-update branch
    obj <- update_m(obj, x, 11)
    expect_equal(dim(obj$S), c(ntemps, p, p))
    expect_true(all(is.finite(obj$S)))

    cand <- gen_cand(obj, x, 12)
    expect_equal(dim(cand), c(ntemps, p))
    expect_true(all(is.finite(cand)))
  })
}


test_that("update_m leaves state untouched before adaptation starts", {
  set.seed(11)
  x <- array(stats::runif(40 * 2 * 2), dim = c(40, 2, 2))
  obj <- make_cov(2, 2, start_adapt = 20)
  before <- obj

  obj <- update_m(obj, x, 5)
  expect_equal(obj$mu, before$mu)
  expect_equal(obj$S, before$S)
})


test_that("update_m recovers the sample mean and covariance at start_adapt_iter", {
  set.seed(12)
  m <- 200; ntemps <- 2; p <- 2
  x <- array(stats::rnorm(m * ntemps * p), dim = c(m, ntemps, p))
  obj <- make_cov(ntemps, p, start_adapt = m)
  obj <- update_m(obj, x, m)

  for (t in 1:ntemps) {
    expect_equal(obj$mu[t, ], colMeans(x[, t, ]), tolerance = 1e-10)
    expect_equal(obj$cov[t, , ], stats::cov(x[, t, ]), tolerance = 1e-10)
  }
})


test_that("update_tau adapts per temperature from its own acceptance count", {
  obj <- make_cov(3, 2, start_adapt = 10, tau = 0)
  # low / high / high acceptance over the last 100 iterations
  obj$count_100 <- c(5, 40, 40)
  out <- update_tau(obj, 100)

  expect_lt(out$tau[1], 0)   # accepted too rarely -> shrink proposals
  expect_gt(out$tau[2], 0)   # accepted too often  -> widen proposals
  expect_gt(out$tau[3], 0)
  expect_equal(out$count_100, rep(0, 3))  # counter resets each window
})


test_that("update_tau only fires on the 100-iteration boundary", {
  obj <- make_cov(2, 1, start_adapt = 10)
  obj$count_100 <- c(0, 50)
  expect_equal(update_tau(obj, 137)$tau, obj$tau)
  expect_equal(update_tau(obj, 137)$count_100, c(0, 50))
})


test_that("gen_cand proposals scale with the proposal covariance", {
  set.seed(13)
  ntemps <- 2; p <- 2; nmcmc <- 5
  x <- array(0.5, dim = c(nmcmc, ntemps, p))

  tight <- make_cov(ntemps, p, start_var = 1e-10)
  wide  <- make_cov(ntemps, p, start_var = 1e-2)

  set.seed(1); dev_tight <- max(abs(gen_cand(tight, x, 2) - 0.5))
  set.seed(1); dev_wide  <- max(abs(gen_cand(wide,  x, 2) - 0.5))

  expect_lt(dev_tight, dev_wide)
})
