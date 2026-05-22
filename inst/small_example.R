\donttest{
  library(impala)
  library(BASS)
  library(mvBayes)
  library(fdasrvf)
  library(bayesplot)

  # generate functions
  f <- function(x) {
    dnorm(seq(0, 1, length.out = 99),
          sin(2 * pi * x[1] ^ 2) / 4 - x[1] / 10 + 0.5,
          0.05) * x[2]
  }

  n = 100
  nt = 99
  p = 3
  x_train = matrix(runif(n * p), n)
  e = rnorm(n * 99)
  y_train = matrix(0, n, nt)
  for (i in 1:n) {
    y_train[i, ] = f(x_train[i, ])
  }

  # generate obs ------------------------------------------------------------
  x_true = c(0.1028, 0.5930)
  ftilde_obs = f(x_true)
  gam_obs = seq(0, 1, length.out = nt)
  vv_obs = gam_to_v(gam_obs)

  tt = seq(0, 1, length.out = nt)
  out = multiple_align_functions(t(y_train), tt, ftilde_obs, 0.01)
  gam_train = out$warping_functions
  vv_train = gam_to_v(gam_train)
  ftilde_train = out$fn
  qtilde_train = out$qn

  # fit emulator ------------------------------------------------------------
  emu_ftilde = mvBayes(bass, x_train, t(ftilde_train), nBasis=2)

  emu_vv = mvBayes(bass, x_train, t(vv_train), nBasis=2)

  # impala ------------------------------------------------------------------
  input_names = c("theta0", "theta1", "theta2")
  bounds = list()
  bounds[['theta0']] = c(0, 1)
  bounds[['theta1']] = c(0, 1)
  bounds[['theta2']] = c(0, 1)

  setup = CalibSetup(bounds, cf_bounds)

  model_ftilde = ModelmvBayes(emu_ftilde, input_names)
  model_vv = ModelmvBayes(emu_vv, input_names)

  setup = addVecExperiments(setup, t(ftilde_obs), model_ftilde, 0.01, 20, rep(1, nt))
  setup = addVecExperiments(setup, t(vv_obs), model_vv, 0.01, 20, rep(1, nt))
  setup = setTemperatureLadder(setup, 1.05 ^ (0:2))
  setup = setMCMC(setup, 1000, 500, 1, 10)
  out_cal = calibPool(setup)
}

