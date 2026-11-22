# Pooled calibration of material strength models on curves they generated,
# ported from python impala's test_physical_models_vec.py::test_constparams
# and test_physmod_init.py. Bounds are +-2% around the generating values, so
# the sampler should get within 1% of each of them.

bounds_around <- function(params) {
  lapply(params, function(v) sort(c(0.98 * v, 1.02 * v)))
}

calibrate_on_own_curve <- function(model_names, params, consts, constraint, nmcmc = 2000) {
  mm <- do.call(MaterialModel, as.list(model_names))
  mm <- setHistoryVariables(mm, emax = 0.6, edot = 2500e-6, nhist = 100)
  mm <- initializeMaterial(mm, params, consts)
  mm <- initializeState(mm, T = 1000)
  h <- computeStateHistory(mm)
  strains <- h[, "strain", 1]
  yobs <- h[, "stress", 1]

  setup <- CalibSetup(bounds_around(params), constraint)
  model <- ModelMaterialStrength(
    temps = 1000, edots = 2500e-6, consts = consts,
    strain_histories = list(strains),
    flow_stress_model = model_names[["flow_stress_model"]],
    melt_model = model_names[["melt_model"]],
    shear_model = model_names[["shear_modulus_model"]],
    specific_heat_model = model_names[["specific_heat_model"]],
    density_model = model_names[["density_model"]]
  )
  # the model reproduces the curve it was generated from
  expect_equal(as.numeric(evalm(model, params)), yobs, tolerance = 1e-8)

  setup <- addVecExperiments(setup, yobs, model, sd_est = 1, s2_df = 0,
                             s2_ind = rep(1, length(yobs)))
  setup <- setTemperatureLadder(setup, 1.05^(0:19))
  setup <- setMCMC(setup, nmcmc = nmcmc, decor = 100)
  calibPool(setup)
}

expect_hits_truth <- function(out, params, tol) {
  expect_setequal(names(out$theta_native), names(params))
  for (k in names(params)) {
    best <- min(abs(out$theta_native[[k]] / params[[k]] - 1))
    expect_lt(best, tol)
  }
}

ptw_true <- list(theta = 0.1, p = 2, s0 = 0.02, sInf = 0.01, kappa = 0.3,
                 lgamma = -12, y0 = 0.01, yInf = 0.003, y1 = 0.09, y2 = 0.7)

sub_models <- function(flow = "PTW_Yield_Stress", heat = "Constant_Specific_Heat",
                       shear = "Stein_Shear_Modulus", melt = "Constant_Melt_Temperature",
                       dens = "Constant_Density") {
  c(flow_stress_model = flow, specific_heat_model = heat, shear_modulus_model = shear,
    melt_model = melt, density_model = dens)
}


test_that("PTW calibration with some parameters held constant", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(1)
  held <- c("y1", "y2", "kappa", "lgamma")
  params <- ptw_true[setdiff(names(ptw_true), held)]
  consts <- c(list(alpha = 0.2, beta = 0.33, matomic = 45.9, Tmelt0 = 2110,
                   rho0 = 4.419, Cv0 = 0.525e-5, G0 = 0.4, chi = 1, sgB = 6.44e-4),
              ptw_true[held])
  # constraints_ptw finds the held-constant PTW parameters among the constants
  out <- calibrate_on_own_curve(sub_models(), params, consts, constraints_ptw)
  expect_hits_truth(out, params, 5e-3)
})


base_consts <- list(alpha = 0.84, beta = 0.33, matomic = 45.9, chi = 1,
                    rho0 = 4.419, rho_0 = 4.45)

physics_cases <- list(
  list(label = "BGP shear, linear melt/heat/density",
       models = sub_models(heat = "Linear_Specific_Heat", shear = "BGP_PW_Shear_Modulus",
                           melt = "Linear_Melt_Temperature", dens = "Linear_Density"),
       consts = list(G0 = 0.44, gamma_1 = 2.2, gamma_2 = -4.7, q2 = 0.8,
                     c0 = 4.730036e-05, c1 = 1.371e-8, tm0 = -3925.796, tm1 = 1448.2,
                     r0 = 4.426741, r1 = -2.5965e-5)),
  list(label = "linear cold PW shear, quadratic melt/heat/density",
       models = sub_models(heat = "Quadratic_Specific_Heat",
                           shear = "Linear_Cold_PW_Shear_Modulus",
                           melt = "Quadratic_Melt_Temperature", dens = "Quadratic_Density"),
       consts = list(g0 = 0.44, g1 = 0.01, c0 = 4.730036e-05, c1 = 1.371e-8, c2 = 1e-10,
                     tm0 = -3925.796, tm1 = 1448.2, tm2 = 5, r0 = 4.426741,
                     r1 = -2.5965e-5, r2 = -1e-8)),
  list(label = "quadratic cold PW shear, cubic melt/heat/density",
       models = sub_models(heat = "Cubic_Specific_Heat",
                           shear = "Quadratic_Cold_PW_Shear_Modulus",
                           melt = "Cubic_Melt_Temperature", dens = "Cubic_Density"),
       consts = list(g0 = 0.44, g1 = 0.01, g2 = 0.001, c0 = 4.730036e-05, c1 = 1.371e-8,
                     c2 = 1e-10, c3 = 1e-12, tm0 = -3925.796, tm1 = 1448.2, tm2 = 5,
                     tm3 = 5, r0 = 4.426741, r1 = -2.5965e-5, r2 = -1e-8, r3 = -1e-10))
)

for (cs in physics_cases) {
  test_that(paste("PTW calibration:", cs$label), {
    skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
    set.seed(2)
    out <- calibrate_on_own_curve(cs$models, ptw_true, c(base_consts, cs$consts),
                             constraints_ptw)
    expect_hits_truth(out, ptw_true, 1e-2)
  })
}


test_that("Steinberg-Guinan calibration with BGP melt and piecewise specific heat", {
  skip_on_cran()  # long MCMC run; exceeds CRAN check-time budget
  set.seed(3)
  # PTW parameters ride along unused, as in the python test; only y0, n and
  # ymax enter the Stein flow stress (beta is held constant)
  params <- c(ptw_true, list(n = 1.001, ymax = 1e6))
  consts <- c(base_consts,
              list(Tm_0 = 1500, G0 = 0.44, g0 = 0.44, g1 = 0.01, g2 = 0.001,
                   gamma_1 = 2.2, gamma_2 = -4.7, gamma_3 = -4.7, q2 = 0.8, q3 = 0.8,
                   T_t = 500, c0_0 = 4.730036e-05, c0_1 = 4.730036e-05, rho_m = 3.5,
                   r0 = 4.426741, c1_0 = 1.371e-8, c2_0 = 1e-10, c3_0 = 1e-12,
                   c1_1 = 1.371e-8, c2_1 = 1e-10, c3_1 = 1e-12,
                   r1 = -2.5965e-5, r2 = -1e-8, r3 = -1e-10, epsi = 1e2))
  models <- sub_models(flow = "Stein_Flow_Stress", heat = "Piecewise_Cubic_Specific_Heat",
                       shear = "Quadratic_Cold_PW_Shear_Modulus",
                       melt = "BGP_Melt_Temperature", dens = "Cubic_Density")
  out <- calibrate_on_own_curve(models, params, consts, cf_bounds)
  expect_hits_truth(out, params, 1e-2)
})
