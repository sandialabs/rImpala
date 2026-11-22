# Material strength physics models, ported from python impala
# (impala.physics.physical_models_vec and physical_models_functions).
#
# Authors of the original python module: DJ Luscher, Peter Trubey, Devin Francom,
# JeeYeon Plohr, Sky Sjue, Lauren VanDervort, Daniel N Blaschke (LANL).
#
# A material model combines one sub-model of each kind (flow stress, specific
# heat, shear modulus, melt temperature, density, heating). Each sub-model is
# described by an entry of `physics_components`: the calibration parameters it
# uses, the constants it needs, default options, and a `value` function of the
# current parameter values `mp` and material state `st`. All quantities are
# vectors with one entry per simulated history, so many parameter settings and
# experiments run at once.

avogadro = 6.022e23
cbrt_avogadro = avogadro^(1 / 3)

cbrt <- function(x) sign(x) * abs(x)^(1 / 3)

erf <- function(x) 2 * stats::pnorm(x * sqrt(2)) - 1

# a constant broadcast to one value per simulated history
const_vec <- function(value, st) rep(value, length.out = length(st$T))


# ---- standalone physics functions -------------------------------------------

cubic_poly <- function(a0, a1, a2, a3, x) a0 + a1 * x + a2 * x^2 + a3 * x^3

piecewise_cubic_specific_heat <- function(Tt, c00, c01, c10, c11, c20, c21, c30, c31, T) {
  hi = T > Tt
  ifelse(hi, cubic_poly(c01, c11, c21, c31, T), cubic_poly(c00, c10, c20, c30, T))
}

bgp_melt_temperature <- function(Tm0, rhom, gamma1, gamma3, q3, rho) {
  Tm0 * cbrt(rho / rhom) *
    exp(6 * gamma1 * (1 / cbrt(rhom) - 1 / cbrt(rho)) +
          2 * gamma3 / q3 * (rhom^(-q3) - rho^(-q3)))
}

quadratic_cold_pw_shear_modulus <- function(g0, g1, g2, alpha, rho, T, Tmelt) {
  gnow = (g0 + g1 * rho + g2 * rho^2) * (1 - alpha * (T / Tmelt))
  gnow[T >= Tmelt] = 0
  gnow[gnow < 0] = 0
  gnow
}

bgp_pw_shear_modulus <- function(G0, rho_0, gamma_1, gamma_2, q2, alpha, rho, T, Tmelt,
                                 melt_strength = 0) {
  cold_shear = G0 * (rho / rho_0)^(4 / 3) *
    exp(6 * gamma_1 * (1 / cbrt(rho_0) - 1 / cbrt(rho)) +
          2 * gamma_2 / q2 * (rho_0^(-q2) - rho^(-q2)))
  gnow = cold_shear * (1 - alpha * pmin(T / Tmelt, 1))
  if (!is.null(melt_strength)) {
    gnow[T >= Tmelt] = melt_strength
  }
  gnow[gnow < 0] = 0
  gnow
}

stein_shear_modulus <- function(G0, sgB, T, Tmelt) {
  # constant density and pressure, so only the temperature term remains
  gnow = G0 * (1 - sgB * (T - 300))
  gnow[T >= Tmelt] = 0
  gnow[gnow < 0] = 0
  gnow
}

jc_yield_stress <- function(edot, A, B, C, n, m, Tref, edot0, eps, T, Tmelt) {
  th = pmax(0, (T - Tref) / (Tmelt - Tref))
  (A + B * eps^n) * (1 + C * log(edot / edot0)) * (1 - th^m)
}

#' @title Check PTW parameters
#' @description Checks the ordering constraints the Preston-Tonks-Wallace (PTW)
#'   flow stress model places on its parameters.
#'
#' @param s0,sInf,y0,yInf,y1,y2,beta PTW parameters (vectors of equal length)
#'
#' @return a logical vector, `TRUE` where the parameter set is valid
#'
#' @export
#'
ptw_goodparam <- function(s0, sInf, y0, yInf, y1, y2, beta) {
  (sInf < s0) & (yInf < y0) & (y0 < s0) & (yInf < sInf) & (y1 > s0) & (y2 > beta)
}

ptw_yield_stress <- function(p, kappa, s0, sInf, y0, yInf, y1, y2, beta, theta, lgamma,
                             edot, rho0, matomic, shear, eps, T, Tmelt, small = 1e-10) {
  n = max(length(T), length(edot), length(p), length(eps), length(shear))
  r = function(x) rep(x, length.out = n)
  p = r(p); kappa = r(kappa); s0 = r(s0); sInf = r(sInf); y0 = r(y0); yInf = r(yInf)
  y1 = r(y1); y2 = r(y2); beta = r(beta); theta = r(theta); lgamma = r(lgamma)
  edot = r(edot); shear = r(shear); eps = r(eps); T = r(T); Tmelt = r(Tmelt)

  t_hom = T / Tmelt
  # 1/a, where 4/3 pi a^3 is the atomic volume
  ainv = cbrt((4 / 3) * pi * rho0 / matomic)
  # transverse wave velocity up to units
  xfact = sqrt(shear / rho0)
  # PTW characteristic strain rate; only its ratio to edot matters
  xiDot = 0.5 * ainv * xfact * cbrt_avogadro

  log_xid_ed = log(xiDot / edot)
  Erfres = erf(kappa * t_hom * (lgamma + log_xid_ed))

  tau_s = pmax(s0 - (s0 - sInf) * Erfres,
               s0 * exp(beta * (-lgamma - log_xid_ed)))
  ayield = y0 - (y0 - yInf) * Erfres
  byield = y1 * exp(-y2 * (lgamma + log_xid_ed))
  cyield = s0 * exp(-beta * (lgamma + log_xid_ed))
  tau_y = pmax(ayield, pmin(byield, cyield))

  scaled_stress = tau_s
  ind = (p > small) & (abs(tau_s - tau_y) > small)
  if (any(ind)) {
    eArg1 = p[ind] * (tau_s[ind] - tau_y[ind]) / (s0[ind] - tau_y[ind])
    eArg2 = eps[ind] * p[ind] * theta[ind] / (s0[ind] - tau_y[ind]) / (exp(eArg1) - 1)
    theLog = log(1 - (1 - exp(-eArg1)) * exp(-eArg2))
    scaled_stress[ind] = tau_s[ind] + (s0[ind] - tau_y[ind]) * theLog / p[ind]
  }
  ind2 = (p <= small) & (tau_s > tau_y)
  if (any(ind2)) {
    d = tau_s[ind2] - tau_y[ind2]
    scaled_stress[ind2] = tau_s[ind2] - d * exp(-eps[ind2] * theta[ind2] / d)
  }
  scaled_stress * shear * 2
}

stein_flow_stress <- function(y0, beta, n, ymax, G0, epsi, shear, eps, T, Tmelt) {
  nn = max(length(T), length(eps), length(shear), length(y0))
  y0 = rep(y0, length.out = nn)
  ymax = rep(ymax, length.out = nn)
  fnow = rep((1 + beta * (epsi + eps))^n, length.out = nn)
  cond1 = fnow * y0 > ymax
  fnow[cond1] = (ymax / y0)[cond1]
  fnow[rep(T, length.out = nn) > Tmelt] = 0
  y0 * fnow * shear / G0
}


# ---- sub-model registry -----------------------------------------------------

component <- function(kind, params = character(), consts = character(),
                      options = list(), value) {
  list(kind = kind, params = params, consts = consts, options = options, value = value)
}

physics_components <- list(
  # heating from plastic work: Taylor-Quinney coefficient chi above edotcrit_chi
  Taylor_Quinney_Heating = component(
    "heating", params = "chi", options = list(edotcrit_chi = 0),
    value = function(mp, st, edot, dt) {
      st$T + (edot > mp$edotcrit_chi) * mp$chi * st$stress * edot * dt / (st$Cv * st$rho)
    }
  ),

  # specific heat
  Constant_Specific_Heat = component(
    "specific_heat", consts = "Cv0",
    value = function(mp, st) const_vec(mp$Cv0, st)
  ),
  Linear_Specific_Heat = component(
    "specific_heat", consts = c("c0", "c1"),
    value = function(mp, st) cubic_poly(mp$c0, mp$c1, 0, 0, st$T)
  ),
  Quadratic_Specific_Heat = component(
    "specific_heat", consts = c("c0", "c1", "c2"),
    value = function(mp, st) cubic_poly(mp$c0, mp$c1, mp$c2, 0, st$T)
  ),
  Cubic_Specific_Heat = component(
    "specific_heat", consts = c("c0", "c1", "c2", "c3"),
    value = function(mp, st) cubic_poly(mp$c0, mp$c1, mp$c2, mp$c3, st$T)
  ),
  Piecewise_Linear_Specific_Heat = component(
    "specific_heat", consts = c("T_t", "c0_0", "c1_0", "c0_1", "c1_1"),
    value = function(mp, st) {
      piecewise_cubic_specific_heat(mp$T_t, mp$c0_0, mp$c0_1, mp$c1_0, mp$c1_1,
                                    0, 0, 0, 0, st$T)
    }
  ),
  Piecewise_Quadratic_Specific_Heat = component(
    "specific_heat", consts = c("T_t", "c0_0", "c1_0", "c2_0", "c0_1", "c1_1", "c2_1"),
    value = function(mp, st) {
      piecewise_cubic_specific_heat(mp$T_t, mp$c0_0, mp$c0_1, mp$c1_0, mp$c1_1,
                                    mp$c2_0, mp$c2_1, 0, 0, st$T)
    }
  ),
  Piecewise_Cubic_Specific_Heat = component(
    "specific_heat",
    consts = c("T_t", "c0_0", "c1_0", "c2_0", "c3_0", "c0_1", "c1_1", "c2_1", "c3_1"),
    value = function(mp, st) {
      piecewise_cubic_specific_heat(mp$T_t, mp$c0_0, mp$c0_1, mp$c1_0, mp$c1_1,
                                    mp$c2_0, mp$c2_1, mp$c3_0, mp$c3_1, st$T)
    }
  ),

  # density
  Constant_Density = component(
    "density", consts = "rho0",
    value = function(mp, st) const_vec(mp$rho0, st)
  ),
  Linear_Density = component(
    "density", consts = c("r0", "r1"),
    value = function(mp, st) cubic_poly(mp$r0, mp$r1, 0, 0, st$T)
  ),
  Quadratic_Density = component(
    "density", consts = c("r0", "r1", "r2"),
    value = function(mp, st) cubic_poly(mp$r0, mp$r1, mp$r2, 0, st$T)
  ),
  Cubic_Density = component(
    "density", consts = c("r0", "r1", "r2", "r3"),
    value = function(mp, st) cubic_poly(mp$r0, mp$r1, mp$r2, mp$r3, st$T)
  ),

  # melt temperature
  Constant_Melt_Temperature = component(
    "melt", consts = "Tmelt0",
    value = function(mp, st) const_vec(mp$Tmelt0, st)
  ),
  Linear_Melt_Temperature = component(
    "melt", consts = c("tm0", "tm1"),
    value = function(mp, st) cubic_poly(mp$tm0, mp$tm1, 0, 0, st$rho)
  ),
  Quadratic_Melt_Temperature = component(
    "melt", consts = c("tm0", "tm1", "tm2"),
    value = function(mp, st) cubic_poly(mp$tm0, mp$tm1, mp$tm2, 0, st$rho)
  ),
  Cubic_Melt_Temperature = component(
    "melt", consts = c("tm0", "tm1", "tm2", "tm3"),
    value = function(mp, st) cubic_poly(mp$tm0, mp$tm1, mp$tm2, mp$tm3, st$rho)
  ),
  # Burakovsky-Greeff-Preston, doi.org/10.1103/PhysRevB.67.094107
  BGP_Melt_Temperature = component(
    "melt", consts = c("Tm_0", "rho_m", "gamma_1", "gamma_3", "q3"),
    value = function(mp, st) {
      bgp_melt_temperature(mp$Tm_0, mp$rho_m, mp$gamma_1, mp$gamma_3, mp$q3, st$rho)
    }
  ),

  # shear modulus
  Constant_Shear_Modulus = component(
    "shear_modulus", consts = "G0",
    value = function(mp, st) const_vec(mp$G0, st)
  ),
  Linear_Cold_PW_Shear_Modulus = component(
    "shear_modulus", consts = c("g0", "g1", "alpha"),
    value = function(mp, st) {
      quadratic_cold_pw_shear_modulus(mp$g0, mp$g1, 0, mp$alpha, st$rho, st$T, st$Tmelt)
    }
  ),
  Quadratic_Cold_PW_Shear_Modulus = component(
    "shear_modulus", consts = c("g0", "g1", "g2", "alpha"),
    value = function(mp, st) {
      quadratic_cold_pw_shear_modulus(mp$g0, mp$g1, mp$g2, mp$alpha, st$rho, st$T, st$Tmelt)
    }
  ),
  Simple_Shear_Modulus = component(
    "shear_modulus", consts = c("G0", "alpha"),
    value = function(mp, st) mp$G0 * (1 - mp$alpha * (st$T / st$Tmelt))
  ),
  # BGP cold shear with PW temperature dependence,
  # Burakovsky, Greeff, Preston, Phys. Rev. B67 (2003) 094107
  BGP_PW_Shear_Modulus = component(
    "shear_modulus", consts = c("G0", "rho_0", "gamma_1", "gamma_2", "q2", "alpha"),
    options = list(melt_strength = 1e-5),
    value = function(mp, st) {
      bgp_pw_shear_modulus(mp$G0, mp$rho_0, mp$gamma_1, mp$gamma_2, mp$q2, mp$alpha,
                           st$rho, st$T, st$Tmelt, mp$melt_strength)
    }
  ),
  # Steinberg-Guinan at constant density and pressure
  Stein_Shear_Modulus = component(
    "shear_modulus", consts = c("G0", "sgB"),
    value = function(mp, st) stein_shear_modulus(mp$G0, mp$sgB, st$T, st$Tmelt)
  ),

  # flow stress
  Constant_Yield_Stress = component(
    "flow_stress", consts = "yield_stress",
    value = function(mp, st, edot) const_vec(mp$yield_stress, st)
  ),
  # Johnson & Cook, Eng. Fract. Mech. 21 (1985) 31
  JC_Yield_Stress = component(
    "flow_stress", params = c("A", "B", "C", "n", "m"), consts = c("Tref", "edot0"),
    value = function(mp, st, edot) {
      jc_yield_stress(edot, mp$A, mp$B, mp$C, mp$n, mp$m, mp$Tref, mp$edot0,
                      st$strain, st$T, st$Tmelt)
    }
  ),
  # Preston, Tonks & Wallace, J. Appl. Phys. 93 (2003) 211
  PTW_Yield_Stress = component(
    "flow_stress",
    params = c("theta", "p", "s0", "beta", "sInf", "kappa", "lgamma",
               "y0", "yInf", "y1", "y2"),
    consts = c("rho0", "matomic"),
    value = function(mp, st, edot) {
      good = ptw_goodparam(mp$s0, mp$sInf, mp$y0, mp$yInf, mp$y1, mp$y2, mp$beta)
      if (!all(good)) {
        cli::cli_abort("PTW parameters violate the model's constraints.",
                       class = "impala_ptw_constraint_error")
      }
      ptw_yield_stress(mp$p, mp$kappa, mp$s0, mp$sInf, mp$y0, mp$yInf, mp$y1, mp$y2,
                       mp$beta, mp$theta, mp$lgamma, edot, mp$rho0, mp$matomic,
                       st$G, st$strain, st$T, st$Tmelt)
    }
  ),
  # Steinberg-Guinan without the thermally activated part
  Stein_Flow_Stress = component(
    "flow_stress", params = c("y0", "beta", "n", "ymax"), consts = c("G0", "epsi"),
    value = function(mp, st, edot) {
      stein_flow_stress(mp$y0, mp$beta, mp$n, mp$ymax, mp$G0, mp$epsi,
                        st$G, st$strain, st$T, st$Tmelt)
    }
  )
)

physics_kinds <- c(
  flow_stress_model = "flow_stress",
  specific_heat_model = "specific_heat",
  shear_modulus_model = "shear_modulus",
  melt_model = "melt",
  density_model = "density",
  heating_model = "heating"
)


#' @title Available material strength sub-models
#' @description Lists the sub-models that [MaterialModel()] and
#'   [ModelMaterialStrength()] accept for each kind of physics.
#'
#' @param name optional sub-model name; when given, returns its calibration
#'   parameters, required constants and default options instead
#'
#' @return a named list of character vectors (one per argument of
#'   [MaterialModel()]), or, with `name`, a list with elements `params`,
#'   `consts` and `options`
#'
#' @export
#'
#' @examples
#' getoptions_ModelMaterialStrength()
#' getoptions_ModelMaterialStrength("PTW_Yield_Stress")
getoptions_ModelMaterialStrength <- function(name = NULL) {
  if (!is.null(name)) {
    comp = physics_components[[name]]
    if (is.null(comp)) {
      cli::cli_abort("Unknown material sub-model {.val {name}}.")
    }
    return(comp[c("params", "consts", "options")])
  }
  kinds = vapply(physics_components, `[[`, "", "kind")
  lapply(physics_kinds, function(k) names(kinds)[kinds == k])
}


#' @title Material strength model
#' @description Builds a material model from one sub-model of each kind. The
#'   model integrates the material state (temperature, strain, stress, shear
#'   modulus, density) along strain histories, vectorised over histories so
#'   many parameter settings or experiments run at once. See
#'   [getoptions_ModelMaterialStrength()] for the available sub-models.
#'
#' @param flow_stress_model flow stress sub-model name
#' @param specific_heat_model specific heat sub-model name
#' @param shear_modulus_model shear modulus sub-model name
#' @param melt_model melt temperature sub-model name
#' @param density_model density sub-model name
#' @param heating_model plastic work heating sub-model name
#'
#' @return An object of class `MaterialModel`. Use [setHistoryVariables()],
#'   [initializeMaterial()] and [initializeState()] to set it up, then
#'   [computeStateHistory()] to run it.
#'
#' @export
#'
#' @examples
#' mm <- MaterialModel(flow_stress_model = "JC_Yield_Stress")
#' mm <- setHistoryVariables(mm, emax = 0.3, edot = 1e-2, nhist = 100)
#' mm <- initializeMaterial(mm,
#'   parameters = list(A = 9e-4, B = 2.92e-3, C = 0.025, n = 0.31, m = 1.09),
#'   constants = list(rho0 = 8.9375, Cv0 = 0.383e-5, G0 = 0.4578,
#'                    Tmelt0 = 1356, Tref = 298, edot0 = 1e-6, chi = 1))
#' mm <- initializeState(mm, T = 298)
#' res <- computeStateHistory(mm)
#' plot(res[, "strain", 1], res[, "stress", 1], type = "l")
MaterialModel <- function(flow_stress_model = "Constant_Yield_Stress",
                          specific_heat_model = "Constant_Specific_Heat",
                          shear_modulus_model = "Constant_Shear_Modulus",
                          melt_model = "Constant_Melt_Temperature",
                          density_model = "Constant_Density",
                          heating_model = "Taylor_Quinney_Heating") {
  chosen = list(
    flow_stress_model = flow_stress_model,
    specific_heat_model = specific_heat_model,
    shear_modulus_model = shear_modulus_model,
    melt_model = melt_model,
    density_model = density_model,
    heating_model = heating_model
  )
  comps = list()
  for (arg in names(chosen)) {
    nm = chosen[[arg]]
    comp = if (is.character(nm) && length(nm) == 1) physics_components[[nm]] else NULL
    if (is.null(comp) || comp$kind != physics_kinds[[arg]]) {
      cli::cli_abort(c(
        "{.arg {arg}} must name a {physics_kinds[[arg]]} sub-model.",
        i = "Options: {.val {getoptions_ModelMaterialStrength()[[arg]]}}."
      ))
    }
    comps[[physics_kinds[[arg]]]] = comp
  }

  # same order as python: flow stress, specific heat, shear, melt, density, heating
  params = unlist(lapply(comps, `[[`, "params"), use.names = FALSE)
  consts = unique(unlist(lapply(comps, `[[`, "consts"), use.names = FALSE))
  options = do.call(c, unname(lapply(comps, `[[`, "options")))
  if (anyDuplicated(params)) {
    cli::cli_abort("Some parameters are shared between sub-models.")
  }
  if (length(intersect(params, consts))) {
    cli::cli_abort("A name is both a parameter and a constant.")
  }

  out = list(
    components = comps,
    model_names = unlist(chosen),
    params = params,
    consts = if (is.null(consts)) character() else consts,
    options = if (is.null(options)) list() else options,
    values = NULL,
    state = list(T = 300, strain = 0, stress = 0),
    strain_history = NULL
  )
  class(out) = "MaterialModel"
  out
}


#' @title Generate strain histories
#' @description Builds evenly spaced constant-rate strain histories from zero
#'   to `emax` at strain rate `edot`.
#'
#' @param emax maximum strain, a scalar or one value per history
#' @param edot strain rate, one value per history
#' @param nhist number of points in each history
#'
#' @return a list with `times`, `strains` (both histories x `nhist`) and
#'   `strain_rate` (histories x (`nhist` - 1))
#'
#' @export
#'
generateStrainHistory <- function(emax, edot, nhist) {
  nexp = max(length(emax), length(edot))
  emax = rep(emax, length.out = nexp)
  edot = rep(edot, length.out = nexp)
  frac = (seq_len(nhist) - 1) / (nhist - 1)
  strains = outer(emax, frac)
  times = outer(emax / edot, frac)
  strains[, nhist] = emax
  times[, nhist] = emax / edot
  rates = (strains[, -1, drop = FALSE] - strains[, -nhist, drop = FALSE]) /
    (times[, -1, drop = FALSE] - times[, -nhist, drop = FALSE])
  list(times = times, strains = strains, strain_rate = rates)
}


#' @title Set strain history of a material model
#' @description Stores constant-rate strain histories for
#'   [computeStateHistory()]; see [generateStrainHistory()].
#'
#' @param obj a `MaterialModel`
#' @param emax maximum strain, a scalar or one value per history
#' @param edot strain rate, one value per history
#' @param nhist number of points in each history
#'
#' @return the updated `MaterialModel`
#'
#' @export
#'
setHistoryVariables <- function(obj, emax, edot, nhist) {
  obj$emax = emax
  obj$edot = edot
  obj$Nhist = nhist
  obj$strain_history = generateStrainHistory(emax, edot, nhist)
  obj
}


#' @title Set parameters of a material model
#' @description Sets the calibration parameters, constants and options of a
#'   material model. A model parameter missing from `parameters` is held
#'   constant at its value in `constants`. Entries of `parameters` the model
#'   does not use are ignored.
#'
#' @param obj a `MaterialModel`
#' @param parameters named list of parameter values (scalars or one value per
#'   history)
#' @param constants named list of constants
#' @param options optional named list overriding the sub-models' default
#'   options (e.g. `edotcrit_chi`, `melt_strength`)
#'
#' @return the updated `MaterialModel`
#'
#' @export
#'
initializeMaterial <- function(obj, parameters, constants, options = NULL) {
  user_constants = setdiff(obj$params, names(parameters))
  need = union(obj$consts, user_constants)
  missing = setdiff(need, names(constants))
  if (length(missing)) {
    cli::cli_abort("Missing value{?s} for {.val {missing}} in {.arg parameters} or {.arg constants}.")
  }
  vals = c(parameters[setdiff(obj$params, user_constants)], constants[need])
  vals = utils::modifyList(vals, obj$options)
  if (!is.null(options)) {
    vals = utils::modifyList(vals, as.list(options))
  }
  obj$values = lapply(vals, as.numeric)
  obj
}


#' @title Set initial state of a material model
#'
#' @param obj a `MaterialModel`
#' @param T initial temperature, a scalar or one value per history
#' @param stress initial stress
#' @param strain initial strain
#'
#' @return the updated `MaterialModel`
#'
#' @export
#'
initializeState <- function(obj, T = 300, stress = 0, strain = 0) {
  obj$state = list(T = T, strain = strain, stress = stress)
  obj
}


# advance the material state by one step of size dt at strain rate edot
update_material_state <- function(obj, st, edot, dt) {
  mp = obj$values
  comps = obj$components
  st$Cv = comps$specific_heat$value(mp, st)
  st$rho = comps$density$value(mp, st)
  st$T = comps$heating$value(mp, st, edot, dt)
  st$strain = st$strain + edot * dt
  st$Tmelt = comps$melt$value(mp, st)
  st$G = comps$shear_modulus$value(mp, st)
  st$stress = comps$flow_stress$value(mp, st, edot)
  st
}


#' @title Run a material model along its strain histories
#'
#' @param obj an initialised `MaterialModel`
#' @param strain_history optional strain history as returned by
#'   [generateStrainHistory()]; defaults to the one set by
#'   [setHistoryVariables()]
#'
#' @return an array of dimension `nhist` x 6 x histories with the time,
#'   strain, stress, temperature, shear modulus and density along each history
#'
#' @export
#'
computeStateHistory <- function(obj, strain_history = NULL) {
  if (is.null(strain_history)) strain_history = obj$strain_history
  if (is.null(strain_history)) {
    cli::cli_abort("No strain history; call {.fn setHistoryVariables} first.")
  }
  if (is.null(obj$values)) {
    cli::cli_abort("Parameters are not set; call {.fn initializeMaterial} first.")
  }
  times = strain_history$times
  rates = strain_history$strain_rate
  nrep = nrow(strain_history$strains)
  nhist = ncol(strain_history$strains)

  out = array(0, c(nhist, 6, nrep),
              dimnames = list(NULL, c("time", "strain", "stress", "T", "G", "rho"), NULL))
  record = function(st, i) {
    out[i, , ] <<- rbind(times[, i],
                         rep(st$strain, length.out = nrep),
                         rep(st$stress, length.out = nrep),
                         rep(st$T, length.out = nrep),
                         rep(st$G, length.out = nrep),
                         rep(st$rho, length.out = nrep))
  }

  st = update_material_state(obj, obj$state, rates[, 1], 0)
  record(st, 1)
  for (i in seq_len(nhist)[-1]) {
    st = update_material_state(obj, st, rates[, i - 1], times[, i] - times[, i - 1])
    record(st, i)
  }
  out
}


#' @title PTW constraint function
#' @description Constraint function for [CalibSetup()] when calibrating the
#'   PTW flow stress model: checks the PTW parameter ordering and the bounds.
#'   PTW parameters held constant are looked up in `constants`.
#'
#' @param x named list of parameter vectors
#' @param bounds list of bounds for each parameter
#' @param constants optional named list of constants, used for PTW parameters
#'   not in `x`. [calibPool()] and friends pass the model's constants.
#'
#' @return a logical vector, `TRUE` where the parameters are valid
#'
#' @export
#'
constraints_ptw <- function(x, bounds, constants = NULL) {
  y = utils::modifyList(as.list(constants), x)
  missing = setdiff(c("s0", "sInf", "y0", "yInf", "y1", "y2", "beta"), names(y))
  if (length(missing)) {
    cli::cli_abort("PTW parameter{?s} {.val {missing}} {?is/are} in neither {.arg x} nor {.arg constants}.")
  }
  good = ptw_goodparam(y$s0, y$sInf, y$y0, y$yInf, y$y1, y$y2, y$beta)
  for (k in names(bounds)) {
    good = good & (x[[k]] < bounds[[k]][2]) & (x[[k]] > bounds[[k]][1])
  }
  as.logical(good)
}
