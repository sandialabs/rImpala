# Material strength physics models, ported from python impala's
# test_physical_models_vec.py, test_physmod_init.py, test_physmod_constyield.py,
# test_physmod_ptw.py and test_physmod_yjc.py.

run_history <- function(mm, params, consts, emax, edot, nhist, T0, options = NULL) {
  mm <- setHistoryVariables(mm, emax = emax, edot = edot, nhist = nhist)
  mm <- initializeMaterial(mm, params, consts, options)
  mm <- initializeState(mm, T = T0, stress = 0, strain = 0)
  computeStateHistory(mm)
}

# strain/stress curve at 2500/s and 1000 K, as in the python regression tests
strainstress <- function(mm, params, consts) {
  h <- run_history(mm, params, consts, emax = 0.6, edot = 2500e-6, nhist = 100, T0 = 1000)
  unname(cbind(h[, "strain", 1], h[, "stress", 1]))
}

read_baseline <- function(i = NULL) {
  f <- if (is.null(i)) "physics_strainstress_baseline.csv" else
    sprintf("physics_strainstress_baseline_%d.csv", i)
  unname(as.matrix(utils::read.csv(test_path("data", f))))
}

ptw_params <- list(theta = 0.1, p = 2, s0 = 0.02, sInf = 0.01, kappa = 0.3,
                   lgamma = -12, y0 = 0.01, yInf = 0.003, y1 = 0.09, y2 = 0.7)

ptw_stein_consts <- list(alpha = 0.2, beta = 0.33, matomic = 45.9, Tmelt0 = 2110,
                         rho0 = 4.419, Cv0 = 0.525e-5, G0 = 0.4, chi = 1, sgB = 6.44e-4)


test_that("PTW with Steinberg-Guinan shear matches the python baseline", {
  mm <- MaterialModel(flow_stress_model = "PTW_Yield_Stress",
                      shear_modulus_model = "Stein_Shear_Modulus")
  expect_equal(strainstress(mm, ptw_params, ptw_stein_consts), read_baseline(),
               tolerance = 1e-8)

  # holding some PTW parameters constant gives the same curve
  held <- c("y1", "y2", "kappa", "lgamma")
  expect_equal(strainstress(mm, ptw_params[setdiff(names(ptw_params), held)],
                            c(ptw_stein_consts, ptw_params[held])),
               read_baseline(), tolerance = 1e-8)
})


test_that("sub-model combinations match the python baselines", {
  base <- list(alpha = 0.84, beta = 0.33, matomic = 45.9, chi = 1, rho0 = 4.419,
               rho_0 = 4.45)

  # 1: BGP shear, linear melt / specific heat / density
  mm <- MaterialModel("PTW_Yield_Stress", "Linear_Specific_Heat", "BGP_PW_Shear_Modulus",
                      "Linear_Melt_Temperature", "Linear_Density")
  consts <- c(base, list(G0 = 0.44, gamma_1 = 2.2, gamma_2 = -4.7, q2 = 0.8,
                         c0 = 4.730036e-05, c1 = 1.371e-8, tm0 = -3925.796, tm1 = 1448.2,
                         r0 = 4.426741, r1 = -2.5965e-5))
  expect_equal(strainstress(mm, ptw_params, consts), read_baseline(1), tolerance = 1e-8)

  # 2: linear cold PW shear, quadratic melt / specific heat / density
  mm <- MaterialModel("PTW_Yield_Stress", "Quadratic_Specific_Heat",
                      "Linear_Cold_PW_Shear_Modulus", "Quadratic_Melt_Temperature",
                      "Quadratic_Density")
  consts <- c(base, list(g0 = 0.44, g1 = 0.01, c0 = 4.730036e-05, c1 = 1.371e-8,
                         c2 = 1e-10, tm0 = -3925.796, tm1 = 1448.2, tm2 = 5,
                         r0 = 4.426741, r1 = -2.5965e-5, r2 = -1e-8))
  expect_equal(strainstress(mm, ptw_params, consts), read_baseline(2), tolerance = 1e-8)

  # 3: quadratic cold PW shear, cubic melt / specific heat / density
  mm <- MaterialModel("PTW_Yield_Stress", "Cubic_Specific_Heat",
                      "Quadratic_Cold_PW_Shear_Modulus", "Cubic_Melt_Temperature",
                      "Cubic_Density")
  consts <- c(base, list(g0 = 0.44, g1 = 0.01, g2 = 0.001, c0 = 4.730036e-05,
                         c1 = 1.371e-8, c2 = 1e-10, c3 = 1e-12, tm0 = -3925.796,
                         tm1 = 1448.2, tm2 = 5, tm3 = 5, r0 = 4.426741,
                         r1 = -2.5965e-5, r2 = -1e-8, r3 = -1e-10))
  expect_equal(strainstress(mm, ptw_params, consts), read_baseline(3), tolerance = 1e-8)

  # 4: Steinberg-Guinan flow stress, BGP melt, piecewise cubic specific heat.
  # beta is a Stein parameter missing from params, so it is held constant.
  mm <- MaterialModel("Stein_Flow_Stress", "Piecewise_Cubic_Specific_Heat",
                      "Quadratic_Cold_PW_Shear_Modulus", "BGP_Melt_Temperature",
                      "Cubic_Density")
  consts <- c(base, list(Tm_0 = 1500, G0 = 0.44, g0 = 0.44, g1 = 0.01, g2 = 0.001,
                         gamma_1 = 2.2, gamma_2 = -4.7, gamma_3 = -4.7, q2 = 0.8,
                         q3 = 0.8, T_t = 500, c0_0 = 4.730036e-05, c0_1 = 4.730036e-05,
                         tm0 = -3925.796, tm1 = 1448.2, tm2 = 5, tm3 = 5, rho_m = 3.5,
                         r0 = 4.426741, c1_0 = 1.371e-8, c2_0 = 1e-10, c3_0 = 1e-12,
                         c1_1 = 1.371e-8, c2_1 = 1e-10, c3_1 = 1e-12,
                         r1 = -2.5965e-5, r2 = -1e-8, r3 = -1e-10, epsi = 1e2))
  params4 <- c(ptw_params, list(n = 1.001, ymax = 1e6))
  expect_equal(strainstress(mm, params4, consts), read_baseline(4), tolerance = 1e-8)
})


test_that("MaterialModel validates sub-models and inputs", {
  expect_error(MaterialModel(flow_stress_model = "Stein_Shear_Modulus"),
               "flow_stress sub-model")
  expect_error(MaterialModel(melt_model = "nope"), "melt sub-model")

  mm <- MaterialModel("PTW_Yield_Stress")
  expect_error(initializeMaterial(mm, ptw_params, list(chi = 1)), "Missing value")
  expect_error(computeStateHistory(mm), "No strain history")

  opts <- getoptions_ModelMaterialStrength()
  expect_true("PTW_Yield_Stress" %in% opts$flow_stress_model)
  expect_true("Stein_Flow_Stress" %in% opts$flow_stress_model)
  expect_true("BGP_Melt_Temperature" %in% opts$melt_model)
  expect_equal(getoptions_ModelMaterialStrength("Stein_Shear_Modulus")$consts,
               c("G0", "sgB"))
})


test_that("PTW refuses parameters that violate its ordering constraints", {
  mm <- MaterialModel("PTW_Yield_Stress", shear_modulus_model = "Stein_Shear_Modulus")
  bad <- ptw_params
  bad$sInf <- 0.03  # sInf must be below s0
  expect_error(strainstress(mm, bad, ptw_stein_consts),
               class = "impala_ptw_constraint_error")

  x <- list(s0 = c(0.02, 0.02), sInf = c(0.01, 0.03))
  consts <- list(y0 = 0.01, yInf = 0.003, y1 = 0.09, y2 = 0.7, beta = 0.33)
  bounds <- list(s0 = c(0, 1), sInf = c(0, 1))
  expect_equal(constraints_ptw(x, bounds, consts), c(TRUE, FALSE))
  expect_equal(constraints_ptw(list(s0 = 2, sInf = 0.01), bounds, consts), FALSE)
})


test_that("histories are vectorised: two runs at once match two single runs", {
  mm <- MaterialModel("PTW_Yield_Stress", shear_modulus_model = "Stein_Shear_Modulus")
  params2 <- lapply(ptw_params, function(v) c(v, v))
  params2$theta <- c(0.1, 0.12)
  h2 <- run_history(mm, params2, ptw_stein_consts, emax = c(0.6, 0.4),
                    edot = c(2500e-6, 1e-4), nhist = 50, T0 = c(1000, 600))
  for (r in 1:2) {
    p1 <- lapply(params2, `[`, r)
    h1 <- run_history(mm, p1, ptw_stein_consts, emax = c(0.6, 0.4)[r],
                      edot = c(2500e-6, 1e-4)[r], nhist = 50, T0 = c(1000, 600)[r])
    expect_equal(h2[, , r], h1[, , 1])
  }
})


# ---- constant yield stress (test_physmod_constyield.py) ---------------------

const_y_consts <- list(yield_stress = 0.1, rho0 = 8.96, Cv0 = 0.383e-5, G0 = 0.4578,
                       Tmelt0 = 1356, chi = 1)

test_that("constant yield: isothermal below edotcrit_chi, adiabatic above", {
  mm <- MaterialModel()
  iso <- run_history(mm, list(), const_y_consts, 1, 1e-7, 1000, 298,
                     options = list(edotcrit_chi = 1e-6))
  expect_equal(iso[1, "T", 1], iso[1000, "T", 1])

  adi <- run_history(mm, list(), const_y_consts, 1, 1, 1000, 298)
  expect_false(adi[1, "T", 1] == adi[1000, "T", 1])
  expect_true(all(adi[, "stress", 1] == const_y_consts$yield_stress))
})


# ---- PTW against FLAG (test_physmod_ptw.py) --------------------------------

ofhc_ptw <- list(theta = 0.025, p = 2, s0 = 0.0085, sInf = 0.00055, kappa = 0.11,
                 lgamma = log(1e-5), y0 = 0.0001, yInf = 0.00009999, y1 = 0.094,
                 y2 = 0.575)
ofhc_base <- list(beta = 0.25, matomic = 63.546, chi = 1, Cv0 = 0.383e-5,
                  rho0 = 8.9375, Tmelt0 = 1625)

ptw_variants <- list(
  list(label = "constant shear", shear = "Constant_Shear_Modulus",
       consts = list(G0 = 0.4578),
       flag = data.frame(emax = c(0.244376, 0.382749, 0.388339, 0.371143),
                         edot = c(1, 1e-2, 1e-4, 1e-6), T0 = 298,
                         stress = c(521.416, 458.198, 418.512, 385.751),
                         tol = c(0.01, 0.09, 0.11, 0.09),
                         opts = c(FALSE, FALSE, FALSE, TRUE))),
  list(label = "simple shear", shear = "Simple_Shear_Modulus",
       consts = list(G0 = 0.50889, alpha = 0.21),
       flag = data.frame(emax = c(0.244463, 0.383542, 0.388302, 0.385416),
                         edot = c(1, 1e-2, 1e-4, 1e-6), T0 = c(620, 298, 298, 298),
                         stress = c(485.301, 439.790, 401.108, 375.227),
                         tol = c(0.06, 0.02, 0.01, 0.02), opts = FALSE)),
  list(label = "Steinberg-Guinan shear", shear = "Stein_Shear_Modulus",
       consts = list(G0 = 0.4578, sgB = 3.8e-4),
       flag = data.frame(emax = c(0.244632, 0.384834, 0.388306, 0.388966),
                         edot = c(1, 1e-2, 1e-4, 1e-6), T0 = c(617, 298, 298, 298),
                         stress = c(427.011, 410.483, 373.629, 351.103),
                         tol = c(0.04, 0.02, 0.01, 0.02), opts = FALSE))
)

rel_diff <- function(a, b) abs(a - b) / mean(c(a, b))

for (v in ptw_variants) {
  test_that(paste("PTW,", v$label, ": isothermal at low rate, adiabatic at high"), {
    mm <- MaterialModel("PTW_Yield_Stress", shear_modulus_model = v$shear)
    consts <- c(ofhc_base, v$consts)
    iso <- run_history(mm, ofhc_ptw, consts, 1, 1e-7, 1000, 298,
                       options = list(edotcrit_chi = 1e-6))
    expect_equal(iso[1, "T", 1], iso[1000, "T", 1])
    adi <- run_history(mm, ofhc_ptw, consts, 1, 1e-2, 1000, 298)
    expect_false(adi[1, "T", 1] == adi[1000, "T", 1])
  })

  test_that(paste("PTW,", v$label, ": flow stress agrees with FLAG"), {
    mm <- MaterialModel("PTW_Yield_Stress", shear_modulus_model = v$shear)
    consts <- c(ofhc_base, v$consts)
    for (k in seq_len(nrow(v$flag))) {
      f <- v$flag[k, ]
      opts <- if (f$opts) list(edotcrit_chi = 1e-6) else NULL
      h <- run_history(mm, ofhc_ptw, consts, f$emax, f$edot, 1000, f$T0, opts)
      expect_lt(rel_diff(h[1000, "stress", 1], f$stress * 1e-5), f$tol)
    }
  })
}


# ---- Johnson-Cook against FLAG (test_physmod_yjc.py) -----------------------

jc_params <- list(A = 0.00090, B = 0.00292, C = 0.0250, n = 0.31, m = 1.09)
jc_consts <- list(rho0 = 8.9375, Cv0 = 0.383e-5, G0 = 0.4578, Tmelt0 = 1356,
                  Tref = 298, edot0 = 1e-6, chi = 1)

test_that("Johnson-Cook: isothermal at low rate, adiabatic at high", {
  mm <- MaterialModel("JC_Yield_Stress")
  iso <- run_history(mm, jc_params, jc_consts, 1, 1e-7, 1000, 298,
                     options = list(edotcrit_chi = 1e-6))
  expect_equal(iso[1, "T", 1], iso[1000, "T", 1])
  adi <- run_history(mm, jc_params, jc_consts, 1, 1e-2, 1000, 298)
  expect_false(adi[1, "T", 1] == adi[1000, "T", 1])
})

test_that("Johnson-Cook flow stress agrees with FLAG", {
  mm <- MaterialModel("JC_Yield_Stress")
  flag <- data.frame(emax = c(0.28179, 0.28895, 0.28928), edot = c(1e-2, 1e-4, 1e-7),
                     stress = c(343.37, 317.06, 284.89), tol = c(0.02, 0.01, 0.05),
                     opts = c(FALSE, FALSE, TRUE))
  for (k in seq_len(nrow(flag))) {
    f <- flag[k, ]
    opts <- if (f$opts) list(edotcrit_chi = 1e-6) else NULL
    h <- run_history(mm, jc_params, jc_consts, f$emax, f$edot, 1000, 298, opts)
    expect_lt(rel_diff(h[1000, "stress", 1], f$stress * 1e-5), f$tol)
  }
})
