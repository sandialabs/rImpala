# Unit tests for the internal helpers whose broadcasting/shape handling the
# calibration relies on. Several of these encode bugs that were silent: the code
# ran and produced numbers, but the numbers were wrong.

test_that("cov_3d_pcm matches a direct per-temperature covariance", {
  cov_3d_pcm <- getFromNamespace("cov_3d_pcm", "impala")

  set.seed(5)
  N <- 50; ntemps <- 3; p <- 2
  arr <- array(stats::rnorm(N * ntemps * p), dim = c(N, ntemps, p))
  mu <- matrix(colMeans(arr), ntemps, p)

  got <- cov_3d_pcm(arr, mu)
  expect_equal(dim(got), c(ntemps, p, p))
  for (t in 1:ntemps) {
    expect_equal(got[t, , ], stats::cov(arr[, t, ]), tolerance = 1e-10)
  }
})


test_that("cov_3d_pcm handles a 1x1 mean without collapsing", {
  # replicate(N, mean, simplify = "array") returns a bare vector when mean is
  # 1 x 1, which used to make aperm() error out for ntemps == p == 1.
  cov_3d_pcm <- getFromNamespace("cov_3d_pcm", "impala")

  set.seed(6)
  arr <- array(stats::rnorm(20), dim = c(20, 1, 1))
  mu <- matrix(mean(arr), 1, 1)
  got <- cov_3d_pcm(arr, mu)

  expect_equal(dim(got), c(1, 1, 1))
  expect_equal(as.numeric(got), stats::var(as.numeric(arr)), tolerance = 1e-10)
})


test_that("s2_kern_sum broadcasts hyperparameters across error groups", {
  s2_kern_sum <- getFromNamespace("s2_kern_sum", "impala")
  ldig_kern <- getFromNamespace("ldig_kern", "impala")

  ntemps <- 3; ns2 <- 2
  ls2 <- matrix(c(1, 2, 3, 10, 20, 30), ntemps, ns2)
  a <- c(0.5, 5); b <- c(1, 100)

  got <- s2_kern_sum(ldig_kern, ls2, a, b)
  # reference: sum over groups, one temperature at a time
  want <- vapply(1:ntemps,
                 function(t) sum(ldig_kern(exp(ls2[t, ]), a, b)),
                 numeric(1))

  expect_length(got, ntemps)
  expect_equal(got, want, tolerance = 1e-10)
  # a plain kern(ls2, a, b) would recycle a/b down temperatures instead
  expect_false(isTRUE(all.equal(got, rowSums(ldig_kern(exp(ls2), a, b)))))
})


test_that("s2_kern_sum and ls2_rowsum accept a single error group", {
  s2_kern_sum <- getFromNamespace("s2_kern_sum", "impala")
  ls2_rowsum <- getFromNamespace("ls2_rowsum", "impala")
  ldig_kern <- getFromNamespace("ldig_kern", "impala")

  ls2 <- matrix(c(0.1, 0.2, 0.3), 3, 1)
  expect_length(s2_kern_sum(ldig_kern, ls2, 2, 1), 3)
  expect_equal(ls2_rowsum(ls2), c(0.1, 0.2, 0.3), tolerance = 1e-12)
})


test_that("chol_solve returns a genuine inverse and log determinant", {
  chol_solve <- getFromNamespace("chol_solve", "impala")

  set.seed(8)
  x <- crossprod(matrix(stats::rnorm(16), 4)) + diag(4)
  got <- chol_solve(x)

  expect_equal(got$inv, solve(x), tolerance = 1e-8)
  expect_equal(got$ldet, as.numeric(determinant(x, logarithm = TRUE)$modulus),
               tolerance = 1e-8)
  expect_equal(got$inv %*% x, diag(4), tolerance = 1e-8)
})


test_that("normalize and unnormalize round-trip", {
  normalize <- getFromNamespace("normalize", "impala")
  unnormalize <- getFromNamespace("unnormalize", "impala")

  bounds <- rbind(c(-2, 5), c(10, 20), c(0, 1))
  z <- matrix(stats::runif(15), 5, 3)

  x <- unnormalize(z, bounds)
  expect_equal(normalize(x, bounds), z, tolerance = 1e-12)
  # each column respects its own bounds -- a recycling bug would mix them
  for (j in 1:3) {
    expect_true(all(x[, j] >= bounds[j, 1] & x[, j] <= bounds[j, 2]))
  }
})


test_that("cf_bounds flags out-of-bounds parameters", {
  bounds <- list(a = c(0, 1), b = c(0, 1))
  inside <- list(a = c(0.5, 0.2), b = c(0.5, 0.9))
  expect_true(all(cf_bounds(inside, bounds)))

  outside <- list(a = c(0.5, 1.5), b = c(0.5, 0.9))
  expect_equal(cf_bounds(outside, bounds), c(TRUE, FALSE))
})


test_that("tran_unif returns a named list on the native scale", {
  bounds_mat <- rbind(c(0, 10), c(-5, 5))
  th <- matrix(c(0, 1, 0.5, 0.5), 2, 2)
  got <- tran_unif(th, bounds_mat, c("alpha", "beta"))

  expect_named(got, c("alpha", "beta"))
  expect_equal(got$alpha, c(0, 10))
  expect_equal(got$beta, c(0, 0))
})
