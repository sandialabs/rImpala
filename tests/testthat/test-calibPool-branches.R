# Tests for calibPool control-flow branches that the shape/correctness suites do
# not reach: the starting-value constraint retry loop and its give-up stop(),
# and the M-H s2 infinity clamp. All use the linear StubModel.

# Minimal setup builder with a custom constraint function.
stub_setup_constrained <- function(constraint_func, p = 2, ny = 8, ntemps = 1,
                                    s2mode = "gibbs", nmcmc = 60, seed = 31) {
  set.seed(seed)
  bounds <- list()
  for (j in 1:p) bounds[[paste0("t_", j)]] <- c(0, 1)
  A <- matrix(stats::rnorm(ny * p), ny, p)
  yobs <- as.numeric(A %*% rep(0.5, p)) + stats::rnorm(ny, 0, 0.05)

  setup <- CalibSetup(bounds, constraint_func)
  setup <- addVecExperiments(setup, yobs, StubModel(A, s2 = s2mode),
                             sd_est = 0.05, s2_df = 20, s2_ind = rep(1, ny))
  setup <- setMCMC(setup, nmcmc = nmcmc, start_adapt_iter = 30, decor = 1000)
  if (ntemps > 1) {
    setup <- setTemperatureLadder(setup, 1.3^(0:(ntemps - 1)), start_temper = 30)
  }
  setup
}


test_that("calibPool retries until the starting constraint is satisfied", {
  # Reject any starting draw whose t_1 is below 0.5, forcing the retry loop to
  # run before it lands a valid start. checkConstraints receives the native
  # named list and must return one logical per temperature.
  constraint <- function(theta, bounds) theta$t_1 >= 0.5

  setup <- stub_setup_constrained(constraint, ntemps = 2)
  set.seed(55)
  out <- suppressMessages(calibPool(setup))

  # the accepted starting value (cold and hot) must satisfy the constraint on
  # the native scale
  start_native <- tran_unif(matrix(out$theta[1, , ], setup$ntemps, setup$p),
                            setup$bounds_mat, names(setup$bounds))
  expect_true(all(start_native$t_1 >= 0.5))
  expect_true(all(is.finite(out$llik)))
})


test_that("calibPool errors when no start can satisfy the constraint", {
  # A constraint that can never be met exhausts the retry budget and stops.
  never <- function(theta, bounds) rep(FALSE, length(theta[[1]]))
  setup <- stub_setup_constrained(never)
  expect_error(
    suppressMessages(calibPool(setup)),
    "starting value satisfying the constraints"
  )
})


test_that("calibPool tolerates an infinite M-H s2 candidate via the clamp", {
  # A wide ls2 proposal can exponentiate to Inf; calibPool clamps that to 1e100
  # in the M-H branch so lik_cov_inv still receives a finite variance. Use a
  # large start_var_ls2 to make oversized proposals likely, and assert the run
  # completes with finite output.
  setup <- stub_setup_constrained(cf_bounds, s2mode = "mh", nmcmc = 200)
  setup$start_var_ls2 <- 25       # inflate the ls2 proposal scale
  set.seed(77)
  out <- suppressMessages(calibPool(setup))

  expect_true(all(is.finite(out$s2[[1]])))
  expect_true(all(is.finite(out$llik)))
})
