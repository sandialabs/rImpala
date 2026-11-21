# Hierarchical and clustered calibration on a linear stub model where each
# group of observations has its own true parameter.

grouped_setup <- function(theta_true, ny_per = 15, s2mode = "gibbs",
                          shared_s2 = FALSE, nmcmc = 1500, ntemps = 1,
                          seed = 5, sd = 0.03) {
  set.seed(seed)
  ng <- nrow(theta_true)
  p <- ncol(theta_true)
  ny <- ng * ny_per
  A <- matrix(stats::rnorm(ny * p), ny, p)
  group <- rep(seq_len(ng), each = ny_per)
  yobs <- rowSums(A * theta_true[group, ]) + stats::rnorm(ny, 0, sd)

  bounds <- list()
  for (k in seq_len(p)) bounds[[paste0("t_", k)]] <- c(0, 1)
  setup <- CalibSetup(bounds, cf_bounds)
  ns2 <- if (shared_s2) 1 else ng
  s2_ind <- if (shared_s2) rep(1, ny) else group
  setup <- addVecExperiments(setup, yobs, StubModel(A, s2 = s2mode),
                             sd_est = rep(sd, ns2), s2_df = rep(20, ns2),
                             s2_ind = s2_ind, theta_ind = group)
  setup <- setMCMC(setup, nmcmc = nmcmc, start_adapt_iter = 200, decor = 50,
                   start_var_theta = 1e-4, start_var_ls2 = 1e-3)
  if (ntemps > 1) {
    setup <- setTemperatureLadder(setup, 1.3^(0:(ntemps - 1)), start_temper = 100)
  }
  setup <- setHierPriors(setup, Sigma0_prior_df = p + 3,
                         Sigma0_prior_scale = diag(p) * 0.05)
  setup
}


test_that("theta_ind is 1-based and complete", {
  A <- matrix(1, 4, 2)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  expect_error(addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4),
                                 theta_ind = c(0, 0, 1, 1)), "1-based")
  expect_error(addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4),
                                 theta_ind = c(1, 1, 3, 3)), "every index")
  s <- addVecExperiments(setup, rep(1, 4), StubModel(A), 0.1, 2, rep(1, 4),
                         theta_ind = c(1, 1, 2, 2))
  expect_equal(s$ntheta, 2)
  expect_equal(s$models[[1]]$nexp, 2)
})


test_that("setHierPriors validates and calibHier requires it", {
  tt <- rbind(c(0.3, 0.6), c(0.5, 0.4))
  setup <- grouped_setup(tt, nmcmc = 5)
  expect_error(setHierPriors(setup, Sigma0_prior_df = 2), "greater than")
  expect_error(setHierPriors(setup, theta0_prior_mean = 0.5), "length 2")
  bare <- setup
  bare$theta0_prior_mean <- NULL
  expect_error(calibHier(bare), "setHierPriors")
  fixed <- fixTheta(setup, "t_1", 0.5)
  expect_error(calibHier(fixed), "not supported")
})


test_that("calibHier recovers group-specific parameters (Gibbs s2)", {
  tt <- rbind(c(0.3, 0.6), c(0.5, 0.4), c(0.7, 0.5))
  setup <- grouped_setup(tt, nmcmc = 2000)
  out <- calibHier(setup)
  keep <- 1001:2000

  expect_equal(dim(out$theta[[1]]), c(2000, 1, 3, 2))
  expect_equal(dim(out$theta0), c(2000, 1, 2))
  expect_equal(dim(out$Sigma0), c(2000, 1, 2, 2))
  for (j in 1:3) {
    expect_equal(colMeans(out$theta[[1]][keep, 1, j, ]), tt[j, ], tolerance = 0.03)
  }
  # theta0 sits between the groups
  expect_true(all(abs(colMeans(out$theta0[keep, 1, ]) - colMeans(tt)) < 0.15))
  expect_equal(nrow(out$theta_parent_native), 2000)
  expect_named(out$theta_native[[1]][[1]], c("t_1", "t_2"))
  # error variance learned near the truth
  expect_true(abs(sqrt(mean(out$s2[[1]][keep, 1, 1])) - 0.03) < 0.015)
})


test_that("calibHier runs with M-H s2 shared across thetas, and tempering", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  tt <- rbind(c(0.3, 0.6), c(0.6, 0.4))
  setup <- grouped_setup(tt, s2mode = "MH", shared_s2 = TRUE, nmcmc = 1500, ntemps = 3)
  out <- calibHier(setup)
  keep <- 801:1500
  expect_equal(dim(out$theta[[1]]), c(1500, 3, 2, 2))
  expect_true(sum(out$count_temper) > 0)
  expect_true(sum(out$count_s2) > 0)
  for (j in 1:2) {
    expect_equal(colMeans(out$theta[[1]][keep, 1, j, ]), tt[j, ], tolerance = 0.05)
  }
  expect_true(all(is.finite(out$llik)))
})


test_that("calibClust groups experiments that share a parameter", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  tt <- rbind(c(0.3, 0.7), c(0.3, 0.7), c(0.7, 0.3), c(0.7, 0.3))
  setup <- grouped_setup(tt, nmcmc = 2000, seed = 8)
  setup <- setClusterPriors(setup, nclustmax = 4, eta_prior_shape = 1, eta_prior_rate = 1)
  out <- calibClust(setup)
  keep <- 1001:2000

  expect_equal(dim(out$theta), c(2000, 1, 4, 2))
  expect_equal(dim(out$delta[[1]]), c(2000, 1, 4))
  d <- out$delta[[1]][keep, 1, ]
  expect_gt(mean(d[, 1] == d[, 2]), 0.8)
  expect_gt(mean(d[, 3] == d[, 4]), 0.8)
  expect_lt(mean(d[, 1] == d[, 3]), 0.2)
  for (j in 1:4) {
    expect_equal(colMeans(out$theta_hist[[1]][keep, 1, j, ]), tt[j, ], tolerance = 0.03)
  }
  expect_true(all(out$eta > 0))
})


test_that("calibClust runs with tempering, M-H s2 and scaled proposals", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  tt <- rbind(c(0.3, 0.7), c(0.7, 0.3), c(0.7, 0.3))
  setup <- grouped_setup(tt, s2mode = "MH", nmcmc = 800, ntemps = 2, seed = 9)
  setup <- setClusterPriors(setup, nclustmax = 3, scale_proposal = TRUE)
  out <- calibClust(setup)
  expect_equal(dim(out$delta[[1]]), c(800, 2, 3))
  expect_true(all(out$delta[[1]] >= 1 & out$delta[[1]] <= 3))
  expect_true(all(is.finite(out$llik)))
  expect_true(sum(out$count_temper) > 0)
})
