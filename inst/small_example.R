library(impala)
library(BASS)
library(mvBayes)

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
f_obs = f(x_true)

tt = seq(0, 1, length.out = nt)

# fit emulator ------------------------------------------------------------
emu = mvBayes(bass, x_train, y_train, nBasis=2)

# impala ------------------------------------------------------------------
input_names = c("theta0", "theta1", "theta2")
bounds = list()
bounds[['theta0']] = c(0, 1)
bounds[['theta1']] = c(0, 1)
bounds[['theta2']] = c(0, 1)

setup = CalibSetup(bounds, cf_bounds)

model = ModelmvBayes(emu, input_names)

setup = addVecExperiments(setup, t(f_obs), model, 0.01, 20, rep(1, nt))
setup = setTemperatureLadder(setup, 1.05 ^ (0:2))
setup = setMCMC(setup, 500, 250, 1, 10)
out_cal = calibPool(setup)

