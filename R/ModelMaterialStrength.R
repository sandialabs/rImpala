# Models with a diagonal Gaussian likelihood: a user-supplied function
# (ModelF) and the material strength simulator (ModelMaterialStrength).

# Diagonal likelihood shared by the models in this file. Weights enter as a
# power likelihood, as in ModelBassPca_func: residuals are scaled by sqrt(wt)
# and each log variance by wt. `inv` is kept as a vector.
diag_llik <- function(yobs, pred, cov, wt = NULL) {
  vec = c(yobs) - c(pred)
  if (is.null(wt)) wt = 1
  -0.5 * cov$ldet - 0.5 * sum(vec * vec * cov$inv * wt)
}

diag_lik_cov_inv <- function(s2vec, wt = NULL) {
  if (is.null(wt)) wt = rep(1, length(s2vec))
  list(inv = 1 / s2vec, ldet = sum(wt * log(s2vec)))
}

diag_discrep_sample <- function(obj, yobs, pred, cov, itemp, wt = NULL) {
  inv = cov$inv
  if (!is.null(wt)) inv = inv * wt
  DtW = t(obj$D * inv)
  S = solve(diag(obj$nd) / obj$discrep_tau + DtW %*% obj$D)
  m = DtW %*% (c(yobs) - c(pred))
  as.numeric(chol_sample(S %*% m, S / itemp))
}

# (nrep x p) matrix of the inputs a model uses, in its order
input_matrix <- function(parmat, input_names) {
  missing = setdiff(input_names, names(parmat))
  if (length(missing)) {
    cli::cli_abort("Calibration parameters lack the model input{?s} {.val {missing}}.")
  }
  do.call(cbind, lapply(input_names, function(v) as.numeric(parmat[[v]])))
}


#' @title Model from a user-supplied function
#'
#' @description Wraps an R function (a simulator or an alternative emulator)
#'   for calibration, with a diagonal Gaussian likelihood.
#'
#' @param f function of one parameter vector, ordered as `input_names`,
#'   returning the prediction of every observation. With `vectorized = TRUE` it
#'   instead takes a matrix with one parameter vector per row and returns a
#'   matrix with one prediction per row.
#' @param input_names names of the calibration parameters `f` takes, in order
#' @param exp_ind experiment indices (default: NULL)
#' @param s2 how to sample the error variance: `"gibbs"` (default) or `"MH"`
#' @param vectorized whether `f` evaluates a matrix of parameters at once
#'
#' @return An object of class `ModelF`
#'
#' @export
#'
#' @examples
#' grid <- seq(0, 1, length.out = 20)
#' model <- ModelF(function(th) th[1] * grid + th[2], c("slope", "intercept"))
ModelF <- function(f, input_names, exp_ind = NULL, s2 = "gibbs", vectorized = FALSE) {
  if (!is.function(f)) {
    cli::cli_abort("{.arg f} must be a function.")
  }
  if (is.null(exp_ind)) exp_ind = 1
  obj = list(
    mod = f,
    input_names = as.character(input_names),
    vectorized = vectorized,
    stochastic = FALSE,
    yobs = NULL,
    meas_error_cor = 1,
    exp_ind = exp_ind,
    nexp = max(exp_ind),
    nd = 0,
    discrep_tau = 1,
    D = NULL,
    s2 = s2,
    constants = NULL,
    supports_wt = TRUE
  )
  class(obj) = "ModelF"
  obj
}


#' @export
evalm.ModelF <- function(obj, parmat, pool = TRUE, nugget = FALSE, ...) {
  x = input_matrix(parmat, obj$input_names)
  if (obj$vectorized) {
    out = obj$mod(x)
    return(matrix(out, nrow(x)))
  }
  rows = lapply(seq_len(nrow(x)), function(r) as.numeric(obj$mod(x[r, ])))
  matrix(unlist(rows), nrow(x), byrow = TRUE)
}

#' @export
llik.ModelF <- function(obj, yobs, pred, cov, wt = NULL, ...) {
  diag_llik(yobs, pred, cov, wt)
}

#' @export
lik_cov_inv.ModelF <- function(obj, s2vec, wt = NULL, inds = NULL, ...) {
  diag_lik_cov_inv(s2vec, wt)
}

#' @export
discrep_sample.ModelF <- function(obj, yobs, pred, cov, itemp, wt = NULL, ...) {
  diag_discrep_sample(obj, yobs, pred, cov, itemp, wt)
}

#' @export
step_m.ModelF <- function(obj, ...) obj


#' @title Material strength model for stress-strain experiments
#'
#' @description Calibrates a [MaterialModel()] against stress-strain curves
#'   from split Hopkinson pressure bar or quasistatic experiments. Each
#'   experiment is simulated at constant strain rate and temperature, and the
#'   simulated stress is interpolated at the measured strains. Observations of
#'   all experiments are stacked in the order of `strain_histories`.
#'
#'   With [calibHier()] or [calibClust()], give each experiment its own
#'   parameter through the `theta_ind` argument of [addVecExperiments()].
#'
#' @param temps temperature of each experiment (Kelvin)
#' @param edots strain rate of each experiment (per microsecond)
#' @param consts named list of constants of the physics sub-models; see
#'   [getoptions_ModelMaterialStrength()] for what each sub-model needs.
#'   Sub-model parameters that are not calibration parameters are taken from
#'   here too.
#' @param strain_histories list of the measured strains of each experiment
#' @param flow_stress_model,melt_model,shear_model,specific_heat_model,density_model
#'   sub-model names, see [getoptions_ModelMaterialStrength()]
#' @param s2 how to sample the error variance: `"gibbs"` (default) or `"MH"`
#' @param options optional named list overriding sub-model options
#' @param nhist number of points in each simulated strain history
#'
#' @return An object of class `ModelMaterialStrength`
#'
#' @export
#'
ModelMaterialStrength <- function(temps,
                                  edots,
                                  consts,
                                  strain_histories,
                                  flow_stress_model,
                                  melt_model,
                                  shear_model,
                                  specific_heat_model,
                                  density_model,
                                  s2 = "gibbs",
                                  options = NULL,
                                  nhist = 100) {
  if (!is.list(strain_histories)) strain_histories = list(strain_histories)
  nexp = length(strain_histories)
  if (length(temps) != nexp || length(edots) != nexp) {
    cli::cli_abort("{.arg temps} and {.arg edots} need one value per strain history ({nexp}).")
  }
  model = MaterialModel(
    flow_stress_model = flow_stress_model,
    specific_heat_model = specific_heat_model,
    shear_modulus_model = shear_model,
    melt_model = melt_model,
    density_model = density_model
  )
  obj = list(
    model = model,
    meas_strain_histories = lapply(strain_histories, as.numeric),
    meas_strain_max = vapply(strain_histories, max, 0),
    constants = consts,
    options = options,
    temps = as.numeric(temps),
    edots = as.numeric(edots),
    nexp_sim = nexp,
    Nhist = nhist,
    stochastic = FALSE,
    yobs = NULL,
    nd = 0,
    discrep_tau = 1,
    D = NULL,
    s2 = s2,
    supports_wt = TRUE
  )
  class(obj) = "ModelMaterialStrength"
  obj
}


#' @export
evalm.ModelMaterialStrength <- function(obj, parmat, pool = TRUE, nugget = FALSE, ...) {
  nrep = length(parmat[[1]])
  nexp = obj$nexp_sim
  # history r = (rep - 1) * nexp + experiment
  parmat_big = lapply(parmat, function(v) rep(as.numeric(v), each = nexp))
  edots = rep(obj$edots, times = nrep)
  temps = rep(obj$temps, times = nrep)
  strain_max = rep(obj$meas_strain_max, times = nrep)
  ntot = nrep * nexp

  mm = setHistoryVariables(obj$model, strain_max, edots, obj$Nhist)
  mm = initializeMaterial(mm, parmat_big, obj$constants, obj$options)
  mm = initializeState(mm, T = temps, stress = rep(0, ntot), strain = rep(0, ntot))
  hist = computeStateHistory(mm)
  # (ntot x nhist)
  sim_strain = t(matrix(hist[, "strain", ], obj$Nhist))
  sim_stress = t(matrix(hist[, "stress", ], obj$Nhist))

  # every rep of an experiment shares its strain grid, so locate the measured
  # strains once per experiment and interpolate all reps together
  pred = vector("list", nexp)
  for (j in seq_len(nexp)) {
    rows = seq(j, ntot, by = nexp)
    grid = sim_strain[rows[1], ]
    x = obj$meas_strain_histories[[j]]
    k = pmin(findInterval(x, grid, rightmost.closed = TRUE, all.inside = TRUE),
             length(grid) - 1)
    w = (x - grid[k]) / (grid[k + 1] - grid[k])
    S = sim_stress[rows, , drop = FALSE]
    lo = S[, k, drop = FALSE]
    hi = S[, k + 1, drop = FALSE]
    pred[[j]] = lo + sweep(hi - lo, 2, w, `*`)
  }
  do.call(cbind, pred)
}

#' @export
llik.ModelMaterialStrength <- function(obj, yobs, pred, cov, wt = NULL, ...) {
  diag_llik(yobs, pred, cov, wt)
}

#' @export
lik_cov_inv.ModelMaterialStrength <- function(obj, s2vec, wt = NULL, inds = NULL, ...) {
  diag_lik_cov_inv(s2vec, wt)
}

#' @export
discrep_sample.ModelMaterialStrength <- function(obj, yobs, pred, cov, itemp, wt = NULL, ...) {
  diag_discrep_sample(obj, yobs, pred, cov, itemp, wt)
}

#' @export
step_m.ModelMaterialStrength <- function(obj, ...) obj
