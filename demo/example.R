library(impala)
library(BASS)
library(mvBayes)
# library(scaledVecchia)
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
x_test = matrix(runif(1000 * p), 1000)
e = rnorm(n * 99)
y_train = matrix(0, n, nt)
for (i in 1:n) {
  y_train[i, ] = f(x_train[i, ])
}
y_test = matrix(0, 1000, nt)
for (i in 1:1000) {
  y_test[i, ] = f(x_test[i, ])
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

matplot(tt, t(y_train), type = "l", col = "gray")
lines(tt, ftilde_obs, col = "black")
legend(
  x = 0,
  y = 8,
  legend = c("simulation", "experiment"),
  col = c("gray", "black"),
  lty = 1
)

matplot(tt, ftilde_train, type = "l", col = "gray")
lines(tt, ftilde_obs, col = "black")
legend(
  x = 0,
  y = 8,
  legend = c("simulation", "experiment"),
  col = c("gray", "black"),
  lty = 1
)

matplot(tt, gam_train, type = "l", col = "gray")
lines(tt, gam_obs, col = "black")
legend(
  x = 0,
  y = 8,
  legend = c("simulation", "experiment"),
  col = c("gray", "black"),
  lty = 1
)

matplot(tt, vv_train, type = "l", col = "gray")
lines(tt, vv_obs, col = "black")
legend(
  x = 0,
  y = 8,
  legend = c("simulation", "experiment"),
  col = c("gray", "black"),
  lty = 1
)

# fit emulator ------------------------------------------------------------
emu_ftilde = mvBayes(bass, x_train, t(ftilde_train), nBasis=2)
plot(emu_ftilde)

emu_vv = mvBayes(bass, x_train, t(vv_train), nBasis=2)
plot(emu_vv)

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
setup = setTemperatureLadder(setup, 1.05 ^ (0:3))
setup = setMCMC(setup, 4000, 2000, 1, 10)
out_cal = calibPool(setup)

# plots -------------------------------------------------------------------
uu = seq(2500, 4000, 2)
theta = tran_unif(out_cal$theta[uu, 1, ], setup$bounds_mat, names(setup$bounds))
expnums = 1

cnt = 1
ftilde_pred_obs = vector(mode = "list", length = setup$nexp)
gam_pred_obs = vector(mode = "list", length = setup$nexp)
for (idx in expnums) {
  time_new = seq(0, 1, length.out = nt)

  fn = setup$models[[1]]$input_names
  parmat_array = matrix(0, length(theta[[fn[1]]]), length(fn))
  for (i in 1:length(fn)) {
    parmat_array[, i] = theta[[fn[i]]]
  }

  ftilde_pred_obs[[idx]] = evalm(setup$models[[cnt]], theta)
  vv_pred_obs = evalm(setup$models[[cnt + 1]], theta)
  cnt = cnt + 2
  gam_pred_obs[[idx]] = v_to_gam(t(vv_pred_obs))

  matplot(
    time_new,
    ftilde_train,
    type = "l",
    col = "gray",
    lty = 1
  )
  matplot(
    time_new,
    t(ftilde_pred_obs[[idx]]),
    type = "l",
    col = "lightblue",
    lty = 1,
    add = TRUE
  )
  lines(time_new, ftilde_obs, col = "black")
  legend(
    x = 0,
    y = 8,
    legend = c("simulation", "prediction", "experiment"),
    col = c("gray", "lightblue", "black"),
    lty = 1
  )

  matplot(
    time_new,
    vv_train,
    type = "l",
    col = "gray",
    lty = 1
  )
  matplot(
    time_new,
    t(vv_pred_obs),
    type = "l",
    col = "lightblue",
    lty = 1,
    add = TRUE
  )
  lines(time_new, vv_obs, col = "black")
  legend(
    x = 0,
    y = 8,
    legend = c("simulation", "prediction", "experiment"),
    col = c("gray", "lightblue", "black"),
    lty = 1
  )

  matplot(
    time_new,
    gam_train,
    type = "l",
    col = "gray",
    lty = 1
  )
  matplot(
    time_new,
    gam_pred_obs[[idx]],
    type = "l",
    col = "lightblue",
    lty = 1,
    add = TRUE
  )
  lines(time_new, gam_obs, col = "black")
  legend(
    x = 0,
    y = 8,
    legend = c("simulation", "prediction", "experiment"),
    col = c("gray", "lightblue", "black"),
    lty = 1
  )
}

samps = data.frame(theta)
mcmc_pairs(samps, off_diag_fun="hex", diag_fun = "dens")

mcmc_hist(samps)

mcmc_trace(samps)

plot(out_cal$llik, type="l", main="Log-Likelihood")

# misaligned prediction
for (idx in expnums){
  obspred = matrix(0, nt, length(uu))
  for (j in 1:length(uu)){
    obspred[,j] = warp_f_gamma(ftilde_pred_obs[[idx]][j,], time_new, invertGamma(gam_pred_obs[[idx]][,j]))
  }

  matplot(
    time_new,
    t(y_train),
    type = "l",
    col = "gray",
    lty = 1
  )
  matplot(
    time_new,
    obspred,
    type = "l",
    col = "lightblue",
    lty = 1,
    add = TRUE
  )
  lines(time_new, ftilde_obs, col = "black")
  legend(
    x = 0,
    y = 8,
    legend = c("simulation", "prediction", "experiment"),
    col = c("gray", "lightblue", "black"),
    lty = 1
  )

}
