library(BASS)
library(fdasrvf)

source("addVecExperiments.R")
source("AMcov_pool.R")
source("calibPool.R")
source("CalibSetup.R")
source("ModelBassPca_func.R")
source("setMCMC.R")
source("setTemperatureLadder.R")
source("UtilityFunctions.R")
source("ClassDefs.R")

# generate functions
f <- function(x) {
  dnorm(seq(0, 1, length.out = 99),
        sin(2 * pi * x[1] ^ 2) / 4 - sqrt(x[1] * x[1]) / 10 + 0.5,
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

# generate obs
x_true = c(0.1028, 0.4930)
ftilde_obs = f(x_true)
gam_obs = seq(0, 1, length.out = nt)
vv_obs = gam_to_v(gam_obs)

tt = seq(0, 1, length.out = nt)
out = multiple_align_functions(t(y_train), tt, ftilde_obs, 0.01)
gam_train = out$warping_functions
vv_train = gam_to_v(gam_train)
ftilde_train = out$fn
qtilde_train = out$qn
ftilde_obs = rowMeans(out$fn)

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


# fit emulator
emu_ftilde = bassPCA(x_train, t(ftilde_train), n.pc = 4)
plot(emu_ftilde)

emu_vv = bassPCA(x_train, t(vv_train), n.pc = 4)
plot(emu_vv)

# impala
input_names = c("theta0", "theta1", "theta2")
bounds = list()
bounds[['theta0']] = c(0, 1)
bounds[['theta1']] = c(0, 1)
bounds[['theta2']] = c(0, 1)

setup = CalibSetup(bounds, cf_bounds)

model_ftilde = ModelBassPca_func(emu_ftilde, input_names)
model_vv = ModelBassPca_func(emu_vv, input_names)

setup = addVecExperiments(setup, t(ftilde_obs), model_ftilde, 0.01, 20, rep(1, nt))
setup = addVecExperiments(setup, t(vv_obs), model_vv, 0.01, 20, rep(1, nt))
setup = setTemperatureLadder(setup, 1.05 ^ (0:3))
setup = setMCMC(setup, 4000, 2000, 1, 10)
out = calibPool(setup)
