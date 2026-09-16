# evalm must return an (ntemps x ny) matrix. R drops dimensions on a
# single-index slice of a 3-d array, so `pred[1, , ]` collapsed to a bare vector
# whenever ntemps == 1 -- the CalibSetup default -- and calibPool then failed
# with "incorrect number of dimensions" while indexing pred_curr[[i]][t, ].
# These are regression tests for that.

test_that("evalm keeps its matrix shape at a single temperature", {
  skip_if_not_installed("BASS")
  # BASS calls parallel::detectCores(), which returns NA in some sandboxed
  # environments and then aborts on `if (n.cores > detectCores())`.
  skip_if(is.na(parallel::detectCores()),
          "parallel::detectCores() is unavailable here")

  set.seed(11)
  nx <- 150; nt <- 10
  tt <- seq(0, 1, length.out = nt)
  X <- matrix(stats::runif(nx * 2), nx, 2)
  fn <- function(x) sin(2 * pi * tt * x[1]) + x[2] * tt
  Y <- t(apply(X, 1, fn))

  bmod <- suppressWarnings(BASS::bassPCA(X, Y, n.pc = 2, nmcmc = 400,
                                         nburn = 300, verbose = FALSE,
                                         n.cores = 1))
  emu <- ModelBassPca_func(bmod, input_names = c("t_1", "t_2"))
  bmat <- rbind(c(0, 1), c(0, 1))

  for (ntemps in c(1, 3)) {
    theta <- matrix(stats::runif(ntemps * 2), ntemps, 2)
    pred <- evalm(emu, tran_unif(theta, bmat, c("t_1", "t_2")), TRUE)
    expect_true(is.matrix(pred))
    expect_equal(dim(pred), c(ntemps, nt))
  }
})


test_that("calibPool runs an emulator untempered and per error group", {
  skip_if_not_installed("BASS")
  skip_if(is.na(parallel::detectCores()),
          "parallel::detectCores() is unavailable here")

  set.seed(11)
  nx <- 150; nt <- 12
  tt <- seq(0, 1, length.out = nt)
  X <- matrix(stats::runif(nx * 2), nx, 2)
  fn <- function(x) sin(2 * pi * tt * x[1]) + x[2] * tt
  Y <- t(apply(X, 1, fn))

  bmod <- suppressWarnings(BASS::bassPCA(X, Y, n.pc = 3, nmcmc = 500,
                                         nburn = 300, verbose = FALSE,
                                         n.cores = 1))
  emu <- ModelBassPca_func(bmod, input_names = c("t_1", "t_2"))

  yobs <- fn(c(0.45, 0.55)) + stats::rnorm(nt, 0, 0.02)
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, yobs, emu, sd_est = c(0.05, 0.05),
                             s2_df = c(2, 2),
                             s2_ind = rep(1:2, each = nt / 2))
  setup <- setMCMC(setup, nmcmc = 600, start_adapt_iter = 200, decor = 100)

  # ntemps defaults to 1, which is exactly the case that used to error
  out <- suppressWarnings(calibPool(setup))
  expect_equal(dim(out$theta), c(600, 1, 2))
  expect_equal(dim(out$s2[[1]]), c(600, 1, 2))
  expect_true(all(is.finite(out$theta)))
  expect_true(all(is.finite(out$s2[[1]])))
  expect_true(all(is.finite(out$llik[-1])))
})
