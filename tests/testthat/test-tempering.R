# Parallel tempering: swap schemes, ladder adaptation, and the consistency of
# each sampler's tempered targets (ports of python impala's test_tempering.py,
# test_clust_tempering.py and test_hier_scale.py).

ladder_setup <- function(tl, scheme = "deo", adapt = FALSE, until = NULL,
                         nmcmc = 2000, nburn = 0) {
  s <- CalibSetup(list(a = c(0, 1)), cf_bounds)
  s <- setMCMC(s, nmcmc = nmcmc, nburn = nburn)
  setTemperatureLadder(s, tl, swap_scheme = scheme, adapt_ladder = adapt,
                       adapt_ladder_until = until)
}


test_that("even/odd pairs and DEO parity alternation", {
  expect_equal(even_odd_pairs(5, 0), rbind(c(1, 3), c(2, 4)))
  expect_equal(even_odd_pairs(5, 1), rbind(c(2, 4), c(3, 5)))
  expect_equal(ncol(even_odd_pairs(1, 0)), 0)
  expect_equal(ncol(even_odd_pairs(2, 1)), 0)

  s <- make_swap_scheduler(ladder_setup(1.1^(0:5)))
  expect_equal(sched_pairs(s)[1, ], c(1, 3, 5))
  expect_equal(sched_pairs(s)[1, ], c(2, 4))
  expect_equal(sched_pairs(s)[1, ], c(1, 3, 5))
})


test_that("random scheme draws disjoint pairs as before", {
  s <- make_swap_scheduler(ladder_setup(1.1^(0:6), "random"), "interleaved")
  set.seed(1)
  sw <- sched_pairs(s)
  set.seed(1)
  old <- matrix(sample(1:7, 6, FALSE), 2)
  expect_equal(sw, old)
  expect_equal(length(unique(as.numeric(sw))), 6)

  s <- make_swap_scheduler(ladder_setup(1.1^(0:6), "random"), "split")
  set.seed(1)
  sw <- sched_pairs(s)
  set.seed(1)
  expect_equal(sw, t(matrix(sample(1:7, 6, FALSE), 3, 2)))
})


test_that("bad tempering options are rejected", {
  s <- CalibSetup(list(a = c(0, 1)), cf_bounds)
  expect_error(setTemperatureLadder(s, 1:3, swap_scheme = "nope"))
  expect_error(setTemperatureLadder(s, 1:3, swap_scheme = "random", adapt_ladder = TRUE),
               "adjacent")
  expect_equal(setTemperatureLadder(s, 1:3)$swap_scheme, "deo")
})


test_that("DEO round-trip rate matches theory and beats SEO", {
  # With every swap accepted DEO carries each replica cold -> hot -> cold in
  # 2(N - 1) rounds, so N replicas complete about rounds / 2 round trips
  # (Syed et al. 2022, Theorem 3 with zero rejection); SEO is diffusive.
  set.seed(0)
  ntemps <- 20
  rounds <- 20000
  trips <- c()
  for (scheme in c("deo", "seo")) {
    s <- make_swap_scheduler(ladder_setup(1.1^(0:(ntemps - 1)), scheme))
    for (r in 1:rounds) {
      sw <- sched_pairs(s)
      sched_record(s, sw, rep(0, ncol(sw)), rep(TRUE, ncol(sw)))
    }
    trips[scheme] <- s$round_trips
  }
  expect_lt(abs(trips[["deo"]] - rounds / 2), 0.05 * rounds / 2)
  expect_lt(trips[["seo"]], trips[["deo"]] / 5)
})


test_that("communication-barrier ladder equalises rejection", {
  # rejection proportional to the gap in inverse temperature -> equal gaps
  itl <- 1 / 1.3^(0:9)
  rej <- 0.5 * -diff(itl)
  new <- communication_barrier_ladder(itl, rej)
  expect_equal(new[1], itl[1])
  expect_equal(new[10], itl[10])
  expect_true(all(diff(new) < 0))
  gaps <- -diff(new)
  expect_equal(gaps, rep(mean(gaps), 9))
})


test_that("ladder adaptation stops at adapt_until", {
  tl <- 1.3^(0:5)
  s <- make_swap_scheduler(ladder_setup(tl, adapt = TRUE, until = 1000))
  before <- s$itl
  changed <- c()
  for (m in 1:1999) {
    sw <- sched_pairs(s)
    # rejection grows towards the cold end
    la <- log(1 - 0.1 * (length(before) - (sw[1, ] - 1)) / length(before))
    sched_record(s, sw, la, rep(FALSE, ncol(sw)))
    if (sched_maybe_adapt(s, m)) changed <- c(changed, m)
  }
  expect_true(length(changed) > 0 && max(changed) <= 1000)
  expect_equal(s$itl[1], before[1])
  expect_equal(s$itl[6], before[6])
  expect_false(isTRUE(all.equal(s$itl, before)))

  # adapt_until defaults to nburn, else nmcmc %/% 2
  expect_equal(make_swap_scheduler(ladder_setup(tl, adapt = TRUE, nmcmc = 900, nburn = 300))$adapt_until, 300)
  expect_equal(make_swap_scheduler(ladder_setup(tl, adapt = TRUE, nmcmc = 900))$adapt_until, 450)
})


test_that("s2 swap kernel matches the s2 update", {
  setup <- list(models = list(list(s2 = "gibbs")), s2_prior_kern = list(ldhc_kern))
  expect_identical(s2_swap_kern(setup, 1), ldig_kern)
  setup$models[[1]]$s2 <- "MH"
  expect_identical(s2_swap_kern(setup, 1), ldhc_kern)
  # with s2_df = 0 the Gibbs update targets p(s2) proportional to 1/s2
  x <- c(0.5, 2)
  expect_equal(ldig_kern(x, 0, 0), -log(x))
})


test_that("sample_eta preserves its Gamma prior (Geweke check)", {
  # if eta ~ Gamma(a, b) and k | eta is the number of clusters of a CRP(eta)
  # partition of n items, one update must return eta ~ Gamma(a, b) again
  set.seed(0)
  a <- 2; b <- 2; n <- 8; N <- 200000
  eta <- stats::rgamma(N, a, rate = b)
  i <- 1:(n - 1)
  k <- 1 + rowSums(matrix(stats::runif(N * (n - 1)), N) < outer(eta, i, function(e, j) e / (e + j)))
  eta_new <- sample_eta(eta, k, n, a, b)
  expect_lt(abs(mean(eta_new) - a / b), 0.01)
  expect_lt(abs(stats::var(eta_new) - a / b^2), 0.02)
  expect_gt(suppressWarnings(stats::ks.test(eta_new[1:20000], "pgamma", a, b))$p.value, 1e-3)
})


gaussian_pool <- function(scheme, adapt) {
  set.seed(3)
  nobs <- 30; sig <- 0.05
  A <- kronecker(diag(2), matrix(1, nobs, 1))
  y <- as.numeric(A %*% c(0.3, 0.6)) + stats::rnorm(2 * nobs, 0, sig)
  setup <- CalibSetup(list(a = c(0, 1), b = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(A, s2 = "fix"), sd_est = sig,
                             s2_df = 0, s2_ind = rep(1, 2 * nobs))
  setup <- setTemperatureLadder(setup, 1.5^(0:7), start_temper = 200,
                                swap_scheme = scheme, adapt_ladder = adapt)
  setup <- setMCMC(setup, nmcmc = 4000, nburn = 1500, decor = 100,
                   start_adapt_iter = 300)
  list(setup = setup, mean = c(mean(y[1:nobs]), mean(y[-(1:nobs)])),
       sd = sig / sqrt(nobs))
}

for (cfg in list(c("deo", "FALSE"), c("deo", "TRUE"), c("random", "FALSE"))) {
  local({
    scheme <- cfg[1]
    adapt <- as.logical(cfg[2])
    test_that(paste0("calibPool recovers a Gaussian posterior (", scheme,
                     if (adapt) " + adapt", ")"), {
      skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
      g <- gaussian_pool(scheme, adapt)
      out <- suppressWarnings(calibPool(g$setup))
      th <- matrix(out$theta[1501:4000, 1, ], ncol = 2)
      expect_true(all(abs(colMeans(th) - g$mean) < 0.25 * g$sd))
      expect_true(all(abs(apply(th, 2, stats::sd) / g$sd - 1) < 0.25))
      expect_length(out$temperature_ladder, 8)
      expect_equal(out$temperature_ladder[1], 1)
      expect_length(out$swap_rejection, 7)
      expect_gte(out$round_trips, 0)
      if (adapt) {
        expect_false(isTRUE(all.equal(out$temperature_ladder, 1.5^(0:7))))
      } else {
        expect_equal(out$temperature_ladder, 1.5^(0:7))
      }
    })
  })
}


test_that("pooled Gibbs s2 posterior is exact with tempering (s2_df = 0)", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  # the model does not depend on theta, so s2 | y ~ IG(n / 2, SS / 2) exactly;
  # the swap must use the 1/s2 kernel the Gibbs update targets
  set.seed(123)
  n <- 20
  y <- stats::rnorm(n, 0, 0.5)
  a <- n / 2; b <- sum(y^2) / 2
  set.seed(1)
  setup <- CalibSetup(list(a = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, y, StubModel(matrix(0, n, 1), s2 = "gibbs"),
                             sd_est = 1, s2_df = 0, s2_ind = rep(1, n))
  setup <- setTemperatureLadder(setup, 2^(0:3), start_temper = 100)
  setup <- setMCMC(setup, nmcmc = 12000, decor = 100)
  out <- suppressWarnings(calibPool(setup))
  s2 <- out$s2[[1]][2001:12000, 1, 1]
  expect_lt(abs(mean(s2) / (b / (a - 1)) - 1), 0.015)
  expect_lt(abs(stats::median(s2) / (b / stats::qgamma(0.5, a)) - 1), 0.015)
})


flat_hier_setup <- function(scale, nmcmc, ntemps = 1, ntheta = 6, nobs = 5) {
  p <- 2
  ind <- rep(seq_len(ntheta), each = nobs)
  setup <- CalibSetup(list(a = c(0, 1), b = c(0, 1)), cf_bounds)
  setup <- addVecExperiments(setup, rep(0, nobs * ntheta),
                             StubModel(matrix(0, nobs * ntheta, p), s2 = "fix"),
                             sd_est = rep(1, ntheta), s2_df = rep(0, ntheta),
                             s2_ind = ind, theta_ind = ind)
  if (ntemps > 1) {
    setup <- setTemperatureLadder(setup, 1.5^(0:(ntemps - 1)), start_temper = 100)
  }
  setup <- setMCMC(setup, nmcmc = nmcmc, decor = 100)
  setHierPriors(setup,
                theta0_prior_mean = rep(if (scale == "probit") 0 else 0.5, p),
                theta0_prior_cov = diag(p) * 0.2^2,
                Sigma0_prior_df = p + 4,
                Sigma0_prior_scale = diag(p) * 0.3^2,
                scale = scale)
}


test_that("setHierPriors validates scale and calibClust refuses probit", {
  setup <- flat_hier_setup("unit", nmcmc = 10)
  expect_equal(setup$hier_scale, "unit")
  expect_error(setHierPriors(setup, scale = "logit"))
  pr <- setHierPriors(setup, scale = "probit")
  expect_equal(pr$theta0_prior_mean, c(0, 0))
  expect_equal(pr$theta0_prior_cov, diag(2))
  expect_error(calibClust(setClusterPriors(pr)), "not supported")
})


for (nt in c(1, 3)) {
  local({
    ntemps <- nt
    test_that(paste0("probit hierarchy recovers the theta0 prior (ntemps = ", ntemps, ")"), {
      skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
      # on the probit scale theta_j ~ N(theta0, Sigma0) is not truncated, so
      # with a flat likelihood theta0 follows its N(0, 0.2^2 I) prior (on the
      # unit scale the bounds tilt it to sd ~ 0.16)
      set.seed(11)
      setup <- flat_hier_setup("probit", nmcmc = 12000, ntemps = ntemps)
      out <- calibHier(setup)
      th0 <- matrix(out$theta0[2001:12000, 1, ], ncol = 2)
      expect_true(all(th0 > 0 & th0 < 1))  # returned on the 0-1 scale
      z0 <- stats::qnorm(clip_unit(th0))
      batches <- apply(array(z0, c(nrow(z0) / 20, 20, 2)), c(2, 3), mean)
      mcse <- apply(batches, 2, stats::sd) / sqrt(20)
      expect_true(all(abs(colMeans(z0)) < 4 * mcse + 0.01))
      expect_true(all(abs(apply(z0, 2, stats::sd) / 0.2 - 1) < 0.12))
      expect_true(all(out$theta[[1]][2001:12000, 1, , ] > 0 &
                        out$theta[[1]][2001:12000, 1, , ] < 1))
      if (ntemps > 1) expect_length(out$temperature_ladder, ntemps)

      parent <- sampleHierParent(setup, th0, out$Sigma0[2001:12000, 1, , ])
      expect_equal(dim(parent), dim(th0))
      expect_true(all(parent > 0 & parent < 1))
    })
  })
}


test_that("sampleHierParent on the unit scale is the old constrained draw", {
  setup <- flat_hier_setup("unit", nmcmc = 10)
  th0 <- matrix(0.5, 50, 2)
  S0 <- array(0, c(50, 2, 2))
  for (k in 1:50) S0[k, , ] <- diag(2) * 0.01
  set.seed(5)
  a <- sampleHierParent(setup, th0, S0)
  set.seed(5)
  b <- t(vapply(1:50, function(k) {
    as.numeric(rmvn_constrained(th0[k, ], t(chol(S0[k, , ])), setup))
  }, numeric(2)))
  expect_identical(a, b)
})


test_that("calibClust with tempering recovers the eta prior under a flat likelihood", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  # with a flat likelihood every chain targets the clustering prior, so the
  # cold chain must reproduce the Gamma(2, 2) prior on eta (mean 1)
  set.seed(1)
  setup <- flat_hier_setup("unit", nmcmc = 5000, ntemps = 3, ntheta = 8, nobs = 4)
  setup <- setClusterPriors(setup, nclustmax = 8, eta_prior_shape = 2, eta_prior_rate = 2)
  out <- calibClust(setup)
  expect_lt(abs(mean(out$eta[1001:5000, 1]) - 1), 0.12)
  expect_length(out$swap_rejection, 2)
})
