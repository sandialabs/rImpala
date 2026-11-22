# Theta priors, ported from python impala's test_theta_prior.py. Input
# validation and per-distribution densities are covered in test-setup.R and
# test-utilities-extra.R; these add the remaining unit checks and the
# end-to-end calibPool checks on a ModelF line model.

eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")

unit_setup <- function() {
  CalibSetup(list(t_0 = c(0, 1), t_1 = c(0, 1)))
}

theta_log_prior <- function(setup, theta) {
  eval_theta_priors(tran_unif(theta, setup$bounds_mat, names(setup$bounds)),
                    setup$theta_prior)
}

prior_cases <- list(
  normal = list(mean = 0.5, sd = 1),
  lognormal = list(meanlog = 0, sdlog = 1),
  beta = list(shape1 = 2, shape2 = 2),
  uniform = list(min = 0, max = 1),
  gamma = list(shape = 2, rate = 3),
  cauchy = list(location = 0.5, scale = 1)
)

for (dist in names(prior_cases)) {
  test_that(paste("the", dist, "prior is finite and integrates to at most 1 on [0, 1]"), {
    setup <- addThetaPrior(unit_setup(), dist, prior_cases[[dist]], "t_0")
    x <- seq(0.01, 0.99, length.out = 25)
    lp <- theta_log_prior(setup, cbind(x, 0.5))
    expect_length(lp, 25)
    expect_true(all(is.finite(lp)))
    dens <- exp(lp)
    trap <- sum(diff(x) * (dens[-1] + dens[-25]) / 2)
    expect_lte(trap, 1 + 1e-8)
  })
}


test_that("a scalar joint prior is applied to every row", {
  setup <- addJointThetaPrior(unit_setup(), c("t_0", "t_1"), function(p) -2.5)
  expect_equal(theta_log_prior(setup, matrix(stats::runif(8), 4)), rep(-2.5, 4))
})


test_that("a prior is -Inf outside its support", {
  setup <- addThetaPrior(unit_setup(), "gamma", list(shape = 2, rate = 1), "t_0")
  expect_equal(theta_log_prior(setup, cbind(0, 0.5)), -Inf)
})


test_that("joint and independent priors add", {
  setup <- addThetaPrior(unit_setup(), "normal", list(mean = 0.5, sd = 0.2), "t_0")
  setup <- addJointThetaPrior(setup, c("t_0", "t_1"), function(p) -1)
  theta <- rbind(c(0.4, 0.6), c(0.5, 0.5))
  expect_equal(theta_log_prior(setup, theta),
               stats::dnorm(theta[, 1], 0.5, 0.2, log = TRUE) - 1)
})


# ---- calibPool with a prior ----------------------------------------------

line_setup <- function(yobs, model, nmcmc = 6000) {
  setup <- CalibSetup(list(t_0 = c(0, 1), t_1 = c(0, 1)))
  setup <- addVecExperiments(setup, yobs, model, sd_est = 0.1, s2_df = 0,
                             s2_ind = rep(1, length(yobs)))
  setup <- setTemperatureLadder(setup, 1.05^(0:7))
  setMCMC(setup, nmcmc = nmcmc, decor = 100)
}

line_model <- function(grid) {
  ModelF(function(th) th[1] * grid + th[2], c("t_0", "t_1"))
}


test_that("a tight prior away from the data pulls the posterior toward it", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(0)
  grid <- seq(0, 1, length.out = 40)
  model <- line_model(grid)
  yobs <- 0.8 * grid + 0.2 + stats::rnorm(40, 0, 0.1)

  set.seed(1)
  out_flat <- calibPool(line_setup(yobs, model))
  setup <- addThetaPrior(line_setup(yobs, model), "normal",
                         list(mean = 0.2, sd = 0.02), "t_0")
  set.seed(1)
  out_prior <- calibPool(setup)

  keep <- 3001:6000
  flat_mean <- mean(out_flat$theta[keep, 1, 1])
  prior_mean <- mean(out_prior$theta[keep, 1, 1])
  expect_lt(abs(flat_mean - 0.8), 0.1)
  expect_lt(prior_mean, flat_mean - 0.2)
  expect_lt(abs(prior_mean - 0.2), 0.05)
})


test_that("with an uninformative likelihood the posterior is the prior", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(0)
  grid <- seq(0, 1, length.out = 20)
  model <- ModelF(function(th) 0 * grid + 0 * th[1], c("t_0", "t_1"))
  yobs <- stats::rnorm(20, 0, 0.1)

  setup <- line_setup(yobs, model, nmcmc = 8000)
  setup <- addThetaPrior(setup, "beta", list(shape1 = 3, shape2 = 6), "t_0")
  set.seed(2)
  out <- calibPool(setup)

  keep <- seq(4001, 8000, by = 2)
  draws <- out$theta[keep, 1, 1]
  a <- 3; b <- 6
  expect_lt(abs(mean(draws) - a / (a + b)), 0.05)
  expect_lt(abs(stats::sd(draws) - sqrt(a * b / ((a + b)^2 * (a + b + 1)))), 0.05)
  # t_1 has no prior, so it stays uniform
  expect_lt(abs(mean(out$theta[keep, 1, 2]) - 0.5), 0.05)
})


test_that("an empty prior list leaves calibPool's draws unchanged", {
  set.seed(0)
  grid <- seq(0, 1, length.out = 30)
  model <- line_model(grid)
  yobs <- 0.8 * grid + 0.2 + stats::rnorm(30, 0, 0.1)

  set.seed(3)
  out_a <- calibPool(line_setup(yobs, model, nmcmc = 500))
  setup <- line_setup(yobs, model, nmcmc = 500)
  setup$theta_prior <- list()
  set.seed(3)
  out_b <- calibPool(setup)
  expect_identical(out_a$theta, out_b$theta)
})
