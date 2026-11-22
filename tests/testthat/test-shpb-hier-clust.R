# Hierarchical and clustered PTW calibration of three SHPB experiments,
# ported from python impala's test_hier.py and test_clust.py. The python tests
# compare the best parent draw to a baseline from numpy's random stream, which
# R cannot reproduce; here the best parent draw must fit the data about as well
# as that python baseline does (within a factor of 2 in SSE).

py_hier_best <- c(0.12125296, 0.889518, 0.7980942, 0.00803078, 0.31638529,
                  0.04834673, 0.49717085, 0.00753855, 0.87773007, 0.03336257)
py_clust_best <- c(0.20526145, 0.59566408, 0.78755387, 0.08794196, 0.33280684,
                   0.44106511, 0.50686367, 0.33024073, 0.64569779, 0.11573176)

best_parent_sse <- function(setup, out, use) {
  parent <- as.matrix(out$theta_parent_native[use, names(setup$bounds)])
  min(shpb_sse(setup, normalize(parent, setup$bounds_mat)))
}


test_that("hierarchical PTW calibration of three SHPB experiments", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(0)
  setup <- shpb_setup(nmcmc = 5000)
  out <- calibHier(setup)
  use <- seq(2501, 5000, by = 2)

  expect_true(all(is.finite(out$llik)))
  expect_lt(best_parent_sse(setup, out, use), 2 * shpb_sse(setup, py_hier_best))

  # each experiment's own parameters fit its curve to well within the noise
  # level the prior expects (sd_est = 1e-4)
  for (j in 1:3) {
    obs <- setup$theta_ind[[1]] == j
    th <- matrix(out$theta[[1]][use, 1, j, ], length(use))
    pred <- evalm(setup$models[[1]], tran_unif(th, setup$bounds_mat, names(setup$bounds)))
    rmse <- sqrt(rowMeans(sweep(pred[, obs, drop = FALSE], 2, setup$ys[[1]][obs])^2))
    expect_lt(stats::median(rmse), 1.5e-4)
  }
})


test_that("clustered PTW calibration of three SHPB experiments", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(1)
  setup <- setClusterPriors(shpb_setup(nmcmc = 2000), nclustmax = 3,
                            eta_prior_shape = 2, eta_prior_rate = 0.1)
  out <- calibClust(setup)
  use <- seq(1001, 2000, by = 2)

  expect_true(all(is.finite(out$llik)))
  expect_lt(best_parent_sse(setup, out, use), 2 * shpb_sse(setup, py_clust_best))
})


test_that("the SHPB model predicts each experiment's stacked observations", {
  setup <- shpb_setup(nmcmc = 10)
  th <- rbind(shpb_best_pool, py_hier_best)
  pred <- evalm(setup$models[[1]], tran_unif(th, setup$bounds_mat, names(setup$bounds)))
  expect_equal(dim(pred), c(2, length(setup$ys[[1]])))
  expect_true(all(is.finite(pred)))
  # the pooled best fit is close to the data everywhere
  expect_lt(max(abs(pred[1, ] - setup$ys[[1]])), 1e-3)
  # rows are evaluated independently
  one <- evalm(setup$models[[1]],
               tran_unif(th[2, , drop = FALSE], setup$bounds_mat, names(setup$bounds)))
  expect_equal(one[1, ], pred[2, ])
  # constraints_ptw rejects the bounds-violating and PTW-violating points
  x <- tran_unif(rbind(shpb_best_pool, c(shpb_best_pool[1:2], 0.01, 0.9,
                                         shpb_best_pool[5:10])),
                 setup$bounds_mat, names(setup$bounds))
  expect_equal(constraints_ptw(x, setup$bounds, setup$constants), c(TRUE, FALSE))
  # beta is a constant here, so leaving out the constants is an error
  expect_error(constraints_ptw(x, setup$bounds), "beta")
})
