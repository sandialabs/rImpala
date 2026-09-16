# Every lik_cov_inv method must return a true inverse and matching log
# determinant. Two of them passed the full matrix to chol2inv() instead of its
# Cholesky factor, which returns a plausible-looking but wrong matrix -- so these
# assert against solve()/determinant() rather than a recorded value.

mk_model <- function(cls, ny, npc = 2) {
  set.seed(15)
  basis <- matrix(stats::rnorm(ny * npc), ny, npc)
  obj <- list(
    meas_error_corr = diag(ny),
    trunc_error_var = diag(ny) * 1e-6,
    discrep_cov = diag(ny) * 1e-12,
    basis = basis,
    emu_vars = rep(1e-4, npc),
    npc = npc
  )
  class(obj) <- cls
  obj
}

classes <- c("ModelmvBayes", "ModelmvBayes_GP",
             "ModelmvBayes_elastic", "ModelmvBayes_elastic_GP")

for (cls in classes) {
  test_that(paste("lik_cov_inv returns a true inverse and ldet:", cls), {
    ny <- 6
    obj <- mk_model(cls, ny)
    s2vec <- seq(0.01, 0.06, length.out = ny)

    got <- lik_cov_inv(obj, s2vec)

    # rebuild the same covariance the methods assemble
    Sigma <- getFromNamespace("cor2cov", "impala")(obj$meas_error_corr, sqrt(s2vec))
    mat <- Sigma + obj$trunc_error_var + obj$discrep_cov +
      obj$basis %*% diag(obj$emu_vars, nrow = obj$npc) %*% t(obj$basis)

    expect_named(got, c("inv", "ldet"))
    expect_equal(got$inv, solve(mat), tolerance = 1e-6)
    expect_equal(got$inv %*% mat, diag(ny), tolerance = 1e-6)
    expect_equal(got$ldet,
                 as.numeric(determinant(mat, logarithm = TRUE)$modulus),
                 tolerance = 1e-8)
    expect_true(isSymmetric(unname(round(got$inv, 8))))
  })
}


for (cls in classes) {
  test_that(paste("lik_cov_inv survives a near-singular covariance:", cls), {
    ny <- 5
    obj <- mk_model(cls, ny)
    # a zero measurement variance makes the plain Cholesky fail; the jitter
    # fallback should keep the method usable
    s2vec <- rep(0, ny)
    obj$trunc_error_var <- matrix(0, ny, ny)
    obj$discrep_cov <- matrix(0, ny, ny)

    got <- lik_cov_inv(obj, s2vec)
    expect_true(all(is.finite(got$inv)))
    expect_true(is.finite(got$ldet))
  })
}


test_that("cor2cov rescales a correlation matrix by the given sds", {
  cor2cov <- getFromNamespace("cor2cov", "impala")
  V <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  sd <- c(2, 3)

  got <- cor2cov(V, sd)
  expect_equal(diag(got), sd^2)
  expect_equal(got[1, 2], 0.5 * 2 * 3)
  # correlations are preserved
  expect_equal(got[1, 2] / prod(sd), V[1, 2])
})


test_that("prior kernels have the documented functional form", {
  ldig <- getFromNamespace("ldig_kern", "impala")
  ldhc <- getFromNamespace("ldhc_kern", "impala")

  x <- c(0.5, 2, 10)
  expect_equal(ldig(x, 2, 3), (-2 - 1) * log(x) - 3 / x)
  expect_equal(ldhc(x, NA, NA), -log(x + 1))
  # inverse-gamma penalises large variances more sharply than half-Cauchy
  expect_lt(ldig(100, 2, 3), ldhc(100, NA, NA))
})
