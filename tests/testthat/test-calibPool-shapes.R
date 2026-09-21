# calibPool has to cope with every combination of p, ns2 and ntemps. Because a
# single-index slice of a 3-d array silently drops dimensions in R, the
# degenerate cases (p == 1, ns2 == 1, ntemps == 1 -- the CalibSetup default) are
# the ones that historically broke. Each case below errored before the indexing
# fixes, so these double as regression tests.

cases <- list(
  list(label = "p=1 ns2=1 ntemps=1 gibbs",    p = 1, ny = 12, ns2 = 1, ntemps = 1, mode = "gibbs"),
  list(label = "p=3 ns2=1 ntemps=1 gibbs",    p = 3, ny = 20, ns2 = 1, ntemps = 1, mode = "gibbs"),
  list(label = "p=3 ns2=1 ntemps=4 tempered", p = 3, ny = 20, ns2 = 1, ntemps = 4, mode = "gibbs"),
  list(label = "p=3 ns2=4 ntemps=4 multi-s2", p = 3, ny = 20, ns2 = 4, ntemps = 4, mode = "gibbs"),
  list(label = "p=2 ns2=3 ntemps=3 mh",       p = 2, ny = 18, ns2 = 3, ntemps = 3, mode = "mh"),
  list(label = "p=1 ns2=1 ntemps=1 mh",       p = 1, ny = 12, ns2 = 1, ntemps = 1, mode = "mh"),
  # ns2 == ny: a separate variance for every component of yobs
  list(label = "p=2 ns2=ny=12 ntemps=1 gibbs", p = 2, ny = 12, ns2 = 12, ntemps = 1, mode = "gibbs"),
  list(label = "p=2 ns2=ny=12 ntemps=3 mh",    p = 2, ny = 12, ns2 = 12, ntemps = 3, mode = "mh")
)

for (cs in cases) {
  test_that(paste("calibPool runs and recovers theta:", cs$label), {
    fx <- stub_setup(cs$p, cs$ny, cs$ns2, cs$ntemps, s2mode = cs$mode)
    out <- suppressWarnings(calibPool(fx$setup))

    # shapes
    expect_equal(dim(out$theta), c(fx$nmcmc, cs$ntemps, cs$p))
    expect_equal(dim(out$s2[[1]]), c(fx$nmcmc, cs$ntemps, cs$ns2))

    # theta_native must be the full cold chain, not a leftover proposal
    expect_s3_class(out$theta_native, "data.frame")
    expect_equal(nrow(out$theta_native), fx$nmcmc)
    expect_equal(ncol(out$theta_native), cs$p)
    expect_true(all(is.finite(as.matrix(out$theta_native))))

    # numerics stay well behaved and inside the unit cube
    expect_true(all(is.finite(out$theta)))
    expect_true(all(is.finite(out$s2[[1]])))
    expect_true(all(is.finite(out$llik[-1])))
    expect_true(all(out$theta >= 0 & out$theta <= 1))

    # the chain actually moves, and lands near the truth
    keep <- (fx$nmcmc %/% 2):fx$nmcmc
    post <- matrix(out$theta[keep, 1, ], length(keep), cs$p)
    expect_gt(stats::sd(as.numeric(post)), 0)
    expect_lt(max(abs(colMeans(post) - fx$theta_true)), 0.35)
  })
}


test_that("calibPool accepts a theta prior", {
  fx <- stub_setup(3, 16, 2, 2, prior = TRUE)
  out <- suppressWarnings(calibPool(fx$setup))
  expect_equal(nrow(out$theta_native), fx$nmcmc)
  expect_true(all(is.finite(out$llik[-1])))
  keep <- (fx$nmcmc %/% 2):fx$nmcmc
  post <- matrix(out$theta[keep, 1, ], length(keep), 3)
  expect_lt(max(abs(colMeans(post) - fx$theta_true)), 0.35)
})


test_that("s2_df = 0 selects the half-Cauchy kernel and still samples", {
  fx <- stub_setup(2, 16, 1, 2, s2_df = 0)
  # any(s2_df == 0) picks ldhc_kern, including the single-group case
  expect_identical(fx$setup$s2_prior_kern[[1]],
                   getFromNamespace("ldhc_kern", "impala"))
  out <- suppressWarnings(calibPool(fx$setup))
  expect_true(all(is.finite(out$llik[-1])))
  expect_equal(nrow(out$theta_native), fx$nmcmc)
})


test_that("s2 chains start from sd_est rather than a constant", {
  fx <- stub_setup(2, 16, 2, 3, sd_est = 0.05)
  out <- suppressWarnings(calibPool(fx$setup))
  # iteration 1 of every temperature/group should be sd_est^2, not exp(1)
  expect_equal(as.numeric(out$s2[[1]][1, , ]),
               rep(0.05^2, 3 * 2), tolerance = 1e-10)
})


test_that("sd_lower and sd_upper confine every s2 draw", {
  # Unbounded, this configuration wanders well above 0.09 (measured maxima of
  # 0.32 under M-H and 0.84 under Gibbs), so the bounds genuinely bite here
  # rather than being satisfied by accident.
  for (mode in c("gibbs", "mh")) {
    fx <- stub_setup(2, 20, 2, 3, s2mode = mode, s2_df = 2, sd_est = 0.05,
                     sd_lower = 0.02, sd_upper = 0.09, nmcmc = 800)
    out <- suppressWarnings(calibPool(fx$setup))
    sd_draws <- sqrt(out$s2[[1]])

    expect_true(all(sd_draws >= 0.02 - 1e-12),
                info = paste(mode, "respects sd_lower"))
    expect_true(all(sd_draws <= 0.09 + 1e-12),
                info = paste(mode, "respects sd_upper"))
    # bounding the variance must not break theta recovery
    keep <- (fx$nmcmc %/% 2):fx$nmcmc
    post <- matrix(out$theta[keep, 1, ], length(keep), 2)
    expect_lt(max(abs(colMeans(post) - fx$theta_true)), 0.35)
  }
})


test_that("omitting the sd bounds leaves the sampler untouched", {
  # The defaults must be inert: unset bounds become (0, Inf) and take neither
  # the Gibbs rejection path nor the M-H rejection path.
  for (mode in c("gibbs", "mh")) {
    run <- function(sd_lower, sd_upper) {
      fx <- stub_setup(2, 16, 2, 2, s2mode = mode, nmcmc = 400,
                       sd_lower = sd_lower, sd_upper = sd_upper)
      set.seed(4242)
      suppressWarnings(calibPool(fx$setup))
    }
    default <- run(NULL, NULL)
    explicit <- run(0, Inf)
    expect_identical(default$theta, explicit$theta)
    expect_identical(default$s2, explicit$s2)
  }
})


test_that("multiple experiments are handled independently", {
  set.seed(3)
  p <- 2; ny <- 16
  bounds <- list(t_1 = c(0, 1), t_2 = c(0, 1))
  A1 <- matrix(stats::rnorm(ny * p), ny, p)
  A2 <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- c(0.4, 0.6)

  setup <- CalibSetup(bounds, cf_bounds)
  for (A in list(A1, A2)) {
    y <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, 0.05)
    setup <- addVecExperiments(setup, y, StubModel(A), sd_est = 0.05,
                              s2_df = 20, s2_ind = rep(1, ny))
  }
  setup <- setMCMC(setup, nmcmc = 400, start_adapt_iter = 100, decor = 50)
  setup <- setTemperatureLadder(setup, 1.3^(0:2), start_temper = 150)

  out <- suppressWarnings(calibPool(setup))
  expect_equal(setup$nexp, 2)
  expect_length(out$s2, 2)
  expect_equal(dim(out$count_s2), c(2, 3))
  for (i in 1:2) expect_true(all(is.finite(out$s2[[i]])))
  keep <- 200:400
  expect_lt(max(abs(colMeans(matrix(out$theta[keep, 1, ], length(keep), p)) -
                      theta_true)), 0.35)
})
