# Tests for the setup objects and their argument validation.

test_that("CalibSetup records bounds and sane defaults", {
  bounds <- list(a = c(0, 1), b = c(-1, 3))
  setup <- CalibSetup(bounds, cf_bounds)

  expect_s3_class(setup, "CalibSetup")
  expect_equal(setup$p, 2)
  expect_equal(setup$nexp, 0)
  expect_equal(setup$bounds_mat, rbind(c(0, 1), c(-1, 3)))
  expect_equal(setup$ntemps, 1)
  expect_null(setup$theta_prior)
  # names must survive, they index tran_unif output
  expect_named(setup$bounds, c("a", "b"))
})


test_that("setMCMC and setTemperatureLadder store their arguments", {
  setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)

  setup <- setMCMC(setup, nmcmc = 500, nburn = 100, thin = 2, decor = 25)
  expect_equal(setup$nmcmc, 500)
  expect_equal(setup$nburn, 100)
  expect_equal(setup$thin, 2)
  expect_equal(setup$decor, 25)

  setup <- setTemperatureLadder(setup, c(1, 2, 4), start_temper = 50)
  expect_equal(setup$ntemps, 3)
  expect_equal(setup$itl, 1 / c(1, 2, 4))
  expect_equal(setup$nswap_per, 1)
  expect_equal(setup$start_temper, 50)
})


test_that("addVecExperiments accumulates experiments", {
  set.seed(2)
  ny <- 10
  setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
  A <- matrix(stats::rnorm(ny), ny, 1)

  setup <- addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                             sd_est = 0.1, s2_df = 5, s2_ind = rep(1, ny))
  expect_equal(setup$nexp, 1)
  expect_equal(setup$ns2[[1]], 1)
  expect_equal(setup$ny_s2[[1]], ny)
  expect_equal(setup$ig_a[[1]], 2.5)

  setup <- addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                             sd_est = 0.2, s2_df = 5, s2_ind = rep(1, ny))
  expect_equal(setup$nexp, 2)
  expect_length(setup$ys, 2)
  expect_length(setup$models, 2)
})


test_that("addVecExperiments counts observations per error group", {
  set.seed(4)
  ny <- 12
  setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
  A <- matrix(stats::rnorm(ny), ny, 1)
  s2_ind <- c(rep(1, 5), rep(2, 7))

  setup <- addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                             sd_est = c(0.1, 0.2), s2_df = c(5, 5),
                             s2_ind = s2_ind)
  expect_equal(setup$ns2[[1]], 2)
  expect_equal(setup$ny_s2[[1]], c(5, 7))
})


test_that("s2_df = 0 selects the half-Cauchy kernel for one group too", {
  set.seed(9)
  ny <- 8
  A <- matrix(stats::rnorm(ny), ny, 1)
  ldhc <- getFromNamespace("ldhc_kern", "impala")
  ldig <- getFromNamespace("ldig_kern", "impala")

  mk <- function(s2_df, sd_est, s2_ind) {
    setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
    addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                      sd_est = sd_est, s2_df = s2_df, s2_ind = s2_ind)
  }

  # a single zero is enough -- the old sum(s2_df == 0) > 1 needed two
  expect_identical(mk(0, 0.1, rep(1, ny))$s2_prior_kern[[1]], ldhc)
  expect_identical(mk(c(0, 0), c(0.1, 0.2), rep(1:2, each = 4))$s2_prior_kern[[1]], ldhc)
  expect_identical(mk(5, 0.1, rep(1, ny))$s2_prior_kern[[1]], ldig)
})


test_that("addVecExperiments rejects malformed s2 groupings", {
  set.seed(11)
  ny <- 20
  A <- matrix(stats::rnorm(ny), ny, 1)
  y <- stats::rnorm(ny)
  mk <- function(...) {
    setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
    addVecExperiments(setup, y, StubModel(A), ...)
  }

  # Each of these used to be accepted silently and either dropped observations
  # or corrupted ig_a/ig_b, sometimes only failing much later inside calibPool.
  # The 0-based case is the one Python users hit: it left a group with no data.
  expect_error(mk(sd_est = rep(0.1, 2), s2_df = rep(2, 2),
                  s2_ind = rep(0:1, each = ny / 2)),
               "must be 1-based")
  expect_error(mk(sd_est = rep(0.1, 2), s2_df = rep(2, 2),
                  s2_ind = rep(1:3, length.out = ny)),
               "must lie between 1 and")
  expect_error(mk(sd_est = 0.1, s2_df = 2, s2_ind = rep(1, ny - 5)),
               "every element of")
  expect_error(mk(sd_est = rep(0.1, 3), s2_df = rep(2, 2),
                  s2_ind = rep(1:3, length.out = ny)),
               "one entry per variance group")
  expect_error(mk(sd_est = c(-0.1, 0.1), s2_df = c(2, 2),
                  s2_ind = rep(1:2, each = ny / 2)),
               "strictly positive")
  expect_error(mk(sd_est = 0.1, s2_df = 2, s2_ind = c(1.5, rep(1, ny - 1))),
               "whole numbers")

  # a group with no observations is legal but prior-driven, so warn rather than stop
  expect_warning(mk(sd_est = rep(0.1, 3), s2_df = rep(2, 3),
                    s2_ind = rep(c(1, 3), each = ny / 2)),
                 "no observations")
})


test_that("addVecExperiments supports one variance per component of yobs", {
  set.seed(13)
  ny <- 15
  A <- matrix(stats::rnorm(ny), ny, 1)
  setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)

  setup <- addVecExperiments(setup, stats::rnorm(ny), StubModel(A),
                             sd_est = rep(0.1, ny), s2_df = rep(2, ny),
                             s2_ind = seq_len(ny))
  expect_equal(setup$ns2[[1]], ny)
  # every component is its own group, so each holds exactly one observation
  expect_equal(setup$ny_s2[[1]], rep(1, ny))
  expect_length(setup$ig_a[[1]], ny)
  expect_length(setup$ig_b[[1]], ny)
})


test_that("addVecExperiments validates and defaults the sd bounds", {
  set.seed(17)
  ny <- 12
  A <- matrix(stats::rnorm(ny), ny, 1)
  y <- stats::rnorm(ny)
  mk <- function(...) {
    setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
    addVecExperiments(setup, y, StubModel(A), sd_est = c(0.1, 0.1),
                      s2_df = c(2, 2), s2_ind = rep(1:2, each = ny / 2), ...)
  }

  # unset bounds become the open interval, so sampling is unconstrained
  unbounded <- mk()
  expect_equal(unbounded$sd_lower[[1]], c(0, 0))
  expect_equal(unbounded$sd_upper[[1]], c(Inf, Inf))

  bounded <- mk(sd_lower = c(0.01, 0.02), sd_upper = c(0.5, 0.6))
  expect_equal(bounded$sd_lower[[1]], c(0.01, 0.02))
  expect_equal(bounded$sd_upper[[1]], c(0.5, 0.6))

  expect_error(mk(sd_lower = 0.01), "one entry per variance group")
  expect_error(mk(sd_lower = c(0.5, 0.5), sd_upper = c(0.2, 0.2)),
               "strictly less than")
  expect_error(mk(sd_lower = c(-1, 0), sd_upper = c(1, 1)), "non-negative")
  # the chain starts at sd_est, so a start outside the bounds could never move
  expect_error(mk(sd_lower = c(0.2, 0.2), sd_upper = c(0.5, 0.5)),
               "must lie within")
})


test_that("addThetaPrior validates the parameter name", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)

  # default pname = NULL used to hit `if (logical(0))` -> "argument is of
  # length zero" rather than the intended message
  expect_error(addThetaPrior(setup), "No parameter name given")
  expect_error(addThetaPrior(setup, pname = "nope"), "No parameter name given")

  ok <- addThetaPrior(setup, "normal", list(mean = 0.5, sd = 1), "t_1")
  expect_length(ok$theta_prior, 1)
  expect_equal(ok$theta_prior[[1]]$name, "t_1")
  expect_equal(ok$theta_prior[[1]]$dist, "normal")

  # priors accumulate rather than overwrite
  ok <- addThetaPrior(ok, "beta", list(shape1 = 2, shape2 = 2), "t_2")
  expect_length(ok$theta_prior, 2)
})


test_that("addJointThetaPrior validates names and the density function", {
  setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)

  expect_error(addJointThetaPrior(setup, c("t_1", "bad"), function(p) 0),
               "not in set of input names")
  expect_error(addJointThetaPrior(setup, c("t_1", "t_2"), "not a function"),
               "must be a function")

  ok <- addJointThetaPrior(setup, c("t_1", "t_2"), function(p) 0)
  expect_length(ok$theta_prior, 1)
  expect_equal(ok$theta_prior[[1]]$names, c("t_1", "t_2"))
})


test_that("eval_theta_priors returns zeros when no prior is set", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")
  theta <- list(t_1 = c(0.2, 0.5, 0.8), t_2 = c(0.1, 0.5, 0.9))
  expect_equal(eval_theta_priors(theta, NULL), rep(0, 3))
})


test_that("eval_theta_priors evaluates independent and joint priors per temperature", {
  eval_theta_priors <- getFromNamespace("eval_theta_priors", "impala")
  theta <- list(t_1 = c(0.2, 0.5, 0.8), t_2 = c(0.1, 0.5, 0.9))

  indep <- list(list(name = "t_1", dist = "normal",
                     params = list(mean = 0.5, sd = 0.2)))
  expect_equal(eval_theta_priors(theta, indep),
               stats::dnorm(theta$t_1, 0.5, 0.2, log = TRUE))

  joint <- list(list(names = c("t_1", "t_2"),
                     log_density_fn = function(p) p$t_1 + p$t_2))
  expect_equal(eval_theta_priors(theta, joint), theta$t_1 + theta$t_2)

  # unknown distributions are rejected rather than silently ignored
  bad <- list(list(name = "t_1", dist = "wat", params = list()))
  expect_error(eval_theta_priors(theta, bad), "Unsupported dist")
})
