# Three split Hopkinson pressure bar experiments on Ti-6Al-4V and a PTW
# setup for them, from python impala's test_hier.py / test_clust.py (based on
# examples/ex_shpb_hierarchical.ipynb). Columns: strain, stress (Mbar).

shpb_data <- function() {
  dat0 <- cbind(seq(0.01, 0.1, by = 0.01),
                c(0.0100583, 0.010381, 0.0106112, 0.0108041, 0.0109619, 0.011082,
                  0.0112103, 0.0113348, 0.0114415, 0.0115094))
  dat1 <- cbind(seq(0.01, 0.32, by = 0.01),
                c(0.004462, 0.0045834, 0.0046879, 0.0047613, 0.004824, 0.0048593,
                  0.0049031, 0.0049449, 0.0049828, 0.005015, 0.0050471, 0.0050752,
                  0.0050964, 0.0051326, 0.0051495, 0.0051803, 0.0051816, 0.0052018,
                  0.0052262, 0.0052345, 0.0052651, 0.0052795, 0.005296, 0.0053118,
                  0.0053274, 0.0053431, 0.00536, 0.0053702, 0.0053762, 0.0053926,
                  0.00541, 0.0054108))
  dat2 <- cbind(
    c(0.005682, 0.012309, 0.020317, 0.02867, 0.036673, 0.045022, 0.05337, 0.061718,
      0.070064, 0.078411, 0.087102, 0.095447, 0.103792, 0.112137, 0.120482, 0.129173,
      0.137517, 0.145862, 0.154206, 0.162897, 0.173673, 0.187578, 0.198008, 0.209478,
      0.217127, 0.225816, 0.234158, 0.243543, 0.250844, 0.259188, 0.26753, 0.276219,
      0.28491, 0.293252, 0.301594, 0.310284, 0.318973, 0.327316, 0.336005, 0.346434),
    c(0.0075202, 0.0077424, 0.0078811, 0.0079921, 0.0080891, 0.0081583, 0.0082275,
      0.0082967, 0.008352, 0.0084074, 0.0084349, 0.0084763, 0.0085177, 0.0085591,
      0.0086005, 0.008628, 0.0086555, 0.0086969, 0.0087244, 0.0087519, 0.0087794,
      0.0088206, 0.0088619, 0.0088754, 0.008903, 0.0089166, 0.0089302, 0.0089438,
      0.0089713, 0.0089988, 0.0090124, 0.009026, 0.0090535, 0.0090671, 0.0090807,
      0.0090943, 0.0091079, 0.0091216, 0.0091352, 0.0091626))
  list(dat = list(dat0, dat1, dat2),
       temps = c(573, 1373, 973),       # Kelvin
       edots = c(800, 2500, 2500))      # 1/s
}

shpb_bounds <- list(
  theta = c(0.0001, 0.2), p = c(0.0001, 5), s0 = c(0.0001, 0.05),
  sInf = c(0.0001, 0.05), kappa = c(0.0001, 0.5), lgamma = c(-14, -9),
  y0 = c(0.0001, 0.05), yInf = c(0.0001, 0.01), y1 = c(0.001, 0.1),
  y2 = c(0.33, 1)
)

# pooled calibration result (0-1 scale) used to start the hierarchical runs
shpb_best_pool <- c(0.165374, 0.135012, 0.71787, 0.001748, 0.315716, 0.033456,
                    0.478695, 0.00509, 0.73037, 0.166709)

# one PTW parameter vector per experiment (theta_ind), as in the python tests
shpb_setup <- function(nmcmc, ntemps = 20) {
  d <- shpb_data()
  model <- ModelMaterialStrength(
    temps = d$temps,
    edots = d$edots * 1e-6,
    consts = list(alpha = 0.2, beta = 0.33, matomic = 45.9, Tmelt0 = 2110,
                  rho0 = 4.419, Cv0 = 0.525e-5, G0 = 0.4, chi = 1, sgB = 6.44e-4),
    strain_histories = lapply(d$dat, function(x) x[, 1]),
    flow_stress_model = "PTW_Yield_Stress",
    melt_model = "Constant_Melt_Temperature",
    shear_model = "Simple_Shear_Modulus",
    specific_heat_model = "Constant_Specific_Heat",
    density_model = "Constant_Density"
  )
  yobs <- unlist(lapply(d$dat, function(x) x[, 2]))
  exp_ind <- rep(seq_along(d$dat), vapply(d$dat, nrow, 0L))

  setup <- CalibSetup(shpb_bounds, constraints_ptw)
  setup <- addVecExperiments(setup, yobs, model,
                             sd_est = rep(1e-4, 3), s2_df = rep(15, 3),
                             s2_ind = exp_ind, theta_ind = exp_ind)
  setup <- setTemperatureLadder(setup, 1.05^(0:(ntemps - 1)), start_temper = 2000)
  setup <- setMCMC(setup, nmcmc = nmcmc, decor = 100)
  p <- setup$p
  setHierPriors(setup,
                theta0_prior_mean = shpb_best_pool,
                theta0_prior_cov = diag(p) * 100,
                Sigma0_prior_df = p + 20,
                Sigma0_prior_scale = diag(p) * 0.1^2,
                theta0_start = matrix(shpb_best_pool, ntemps, p, byrow = TRUE))
}

# sum of squared errors of the predictions at each row of a 0-1 scale matrix
shpb_sse <- function(setup, mat) {
  mat <- matrix(mat, ncol = setup$p)
  pred <- evalm(setup$models[[1]], tran_unif(mat, setup$bounds_mat, names(setup$bounds)))
  rowSums(sweep(pred, 2, setup$ys[[1]])^2)
}
