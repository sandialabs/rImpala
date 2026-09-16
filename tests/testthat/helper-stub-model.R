# A minimal linear emulator used to exercise calibPool without depending on
# BASS/mvBayes (both only Suggests). The forward model is y = A %*% theta with a
# diagonal likelihood, so the posterior has a closed form -- see
# test-posterior-correctness.R, which compares against it.

StubModel <- function(A, s2 = "gibbs") {
  ny <- nrow(A)
  obj <- list(
    A = A,
    ny = ny,
    p = ncol(A),
    s2 = s2,
    nd = 0,
    stochastic = FALSE,
    meas_error_corr = diag(ny),
    trunc_error_var = matrix(0, ny, ny),
    discrep_cov = diag(ny) * 1e-12,
    basis = matrix(0, ny, 1),
    emu_vars = 0,
    npc = 1
  )
  class(obj) <- "StubModel"
  obj
}

evalm.StubModel <- function(obj, parmat, pool = TRUE, nugget = FALSE, ...) {
  t(obj$A %*% t(do.call(cbind, parmat)))
}

step_m.StubModel <- function(obj, ...) obj

llik.StubModel <- function(obj, yobs, pred, cov, ...) {
  vec <- c(yobs - pred)
  as.numeric(-0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec))
}

lik_cov_inv.StubModel <- function(obj, s2vec, ...) {
  R <- chol(diag(s2vec, nrow = length(s2vec)) + obj$discrep_cov)
  list(inv = chol2inv(R), ldet = 2 * sum(log(diag(R))))
}

registerS3method("evalm", "StubModel", evalm.StubModel)
registerS3method("step_m", "StubModel", step_m.StubModel)
registerS3method("llik", "StubModel", llik.StubModel)
registerS3method("lik_cov_inv", "StubModel", lik_cov_inv.StubModel)


# Build a CalibSetup around StubModel for a given shape/mode combination.
stub_setup <- function(p, ny, ns2, ntemps, s2mode = "gibbs", nmcmc = 400,
                       s2_df = 20, sd_est = 0.05, prior = FALSE, seed = 42,
                       sd_lower = NULL, sd_upper = NULL) {
  set.seed(seed)
  bounds <- list()
  for (j in 1:p) bounds[[paste0("t_", j)]] <- c(0, 1)

  A <- matrix(stats::rnorm(ny * p), ny, p)
  theta_true <- seq(0.3, 0.7, length.out = p)
  yobs <- as.numeric(A %*% theta_true) + stats::rnorm(ny, 0, sd_est)

  setup <- CalibSetup(bounds, cf_bounds)
  setup <- addVecExperiments(
    setup, yobs, StubModel(A, s2 = s2mode),
    sd_est = rep(sd_est, ns2),
    s2_df  = rep(s2_df, ns2),
    s2_ind = rep(1:ns2, length.out = ny),
    sd_lower = if (is.null(sd_lower)) NULL else rep(sd_lower, ns2),
    sd_upper = if (is.null(sd_upper)) NULL else rep(sd_upper, ns2)
  )
  if (prior) {
    setup <- addThetaPrior(setup, "normal", list(mean = 0.5, sd = 0.3), "t_1")
  }
  setup <- setMCMC(setup, nmcmc = nmcmc, start_adapt_iter = 100, decor = 50)
  if (ntemps > 1) {
    setup <- setTemperatureLadder(setup, 1.3^(0:(ntemps - 1)),
                                  start_temper = 150)
  }
  list(setup = setup, theta_true = theta_true, nmcmc = nmcmc, p = p)
}
