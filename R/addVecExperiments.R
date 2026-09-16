#' @title Add vector experiments
#'
#' @description This method adds vector experiments to calibration object
#'
#' @param obj `CalibSetup` Object
#' @param yobs a vector of the experiment or observation
#' @param model emulator (currently expecting a object of class
#'                         `ModelBassPca_func` or `ModelmvBayes`)
#' @param sd_est vector of initial standard deviation estimates, one per
#'   separately estimated measurement error variance. `length(sd_est)` sets the
#'   number of variance groups.
#' @param s2_df vector of inverse gamma prior degrees of freedom, one per
#'   variance group, so `length(s2_df)` must equal `length(sd_est)`. A value of
#'   zero selects a half-Cauchy prior instead; if any entry is zero the
#'   half-Cauchy is used for every group.
#' @param s2_ind vector of 1-based variance group indices, one per element of
#'   `yobs`, mapping each observation to the variance it is drawn with. Must
#'   satisfy `length(s2_ind) == length(yobs)` and
#'   `max(s2_ind) <= length(sd_est)`.
#' @param sd_lower optional vector of lower bounds on the learned standard
#'   deviations, one per variance group (default: `NULL`, meaning no lower
#'   bound)
#' @param sd_upper optional vector of upper bounds on the learned standard
#'   deviations, one per variance group (default: `NULL`, meaning no upper
#'   bound)
#' @param meas_error_cor measurement error correlation (default: `NULL`)
#' @param theta_ind indices of theta (default: `NULL`)
#' @param D discrepancy basis (matrix of columns of basis, default:
#'           `NaN`)
#' @param discrep_tau discrepancy sampling tau
#'
#' @return An object of class `CalibSetup`
#'
#' @details
#' `sd_est`, `s2_df`, `sd_lower` and `sd_upper` are all indexed by variance
#' group, while `s2_ind` is indexed by observation. A scalar `sd_est` with
#' `s2_ind = rep(1, length(yobs))` gives the usual single measurement error
#' variance for the whole vector.
#'
#' Any grouping of the observations is allowed. Passing
#' `s2_ind = seq_along(yobs)` with `sd_est` of that same length learns a
#' separate standard deviation for every component of `yobs`. Note that each
#' variance is then informed by a single observation, so the prior set through
#' `s2_df` and `sd_est` does most of the work; use a larger `s2_df`, or
#' `sd_lower`/`sd_upper`, when the per component variances need to be
#' constrained.
#'
#' @examples
#' set.seed(1)
#' ny <- 20
#' A <- matrix(stats::rnorm(ny * 2), ny, 2)
#' yobs <- as.numeric(A %*% c(0.45, 0.55)) + stats::rnorm(ny, 0, 0.05)
#'
#' # a minimal linear emulator standing in for a fitted surrogate
#' mod <- structure(list(s2 = "gibbs", nd = 0, stochastic = FALSE,
#'                       discrep_cov = diag(ny) * 1e-12), class = "ExampleModel")
#'
#' setup <- CalibSetup(list(t_1 = c(0, 1), t_2 = c(0, 1)), cf_bounds)
#'
#' # two error groups: the first half of yobs shares one variance, the rest another
#' two_group <- addVecExperiments(setup, yobs, mod,
#'                                sd_est = c(0.05, 0.05),
#'                                s2_df  = c(2, 2),
#'                                s2_ind = rep(1:2, each = ny / 2))
#' two_group$ns2[[1]]
#'
#' # one standard deviation per component of yobs, bounded away from zero
#' per_component <- addVecExperiments(setup, yobs, mod,
#'                                    sd_est   = rep(0.05, ny),
#'                                    s2_df    = rep(10, ny),
#'                                    s2_ind   = seq_along(yobs),
#'                                    sd_lower = rep(0.01, ny),
#'                                    sd_upper = rep(0.50, ny))
#' per_component$ns2[[1]]
#'
#' @export
#'
addVecExperiments <- function(obj,
                              yobs,
                              model,
                              sd_est,
                              s2_df,
                              s2_ind,
                              sd_lower = NULL,
                              sd_upper = NULL,
                              meas_error_cor = NULL,
                              theta_ind = NULL,
                              D = NULL,
                              discrep_tau = 1) {
  N = length(obj$ys)

  # Validate before touching `obj`: a partially updated setup is worse than an
  # error, and several of these mistakes used to corrupt the run silently.
  ns2 = length(sd_est)

  if (!is.numeric(sd_est) || ns2 == 0) {
    cli::cli_abort("{.arg sd_est} must be a non-empty numeric vector, one entry per variance group.")
  }
  if (anyNA(sd_est) || any(!is.finite(sd_est)) || any(sd_est <= 0)) {
    cli::cli_abort("{.arg sd_est} must be finite and strictly positive.")
  }
  if (length(s2_df) != ns2) {
    cli::cli_abort(c(
      "{.arg s2_df} and {.arg sd_est} must have one entry per variance group.",
      i = "{.arg s2_df} has length {length(s2_df)} but {.arg sd_est} has length {ns2}."
    ))
  }
  if (!is.numeric(s2_df) || anyNA(s2_df) || any(!is.finite(s2_df)) || any(s2_df < 0)) {
    cli::cli_abort("{.arg s2_df} must be finite and non-negative.")
  }
  if (!is.numeric(s2_ind) || anyNA(s2_ind)) {
    cli::cli_abort("{.arg s2_ind} must be a numeric vector with no missing values.")
  }
  if (length(s2_ind) != length(yobs)) {
    cli::cli_abort(c(
      "{.arg s2_ind} must give a variance group for every element of {.arg yobs}.",
      i = "{.arg s2_ind} has length {length(s2_ind)} but {.arg yobs} has length {length(yobs)}."
    ))
  }
  if (any(s2_ind != as.integer(s2_ind))) {
    cli::cli_abort("{.arg s2_ind} must contain whole numbers.")
  }
  if (any(s2_ind == 0)) {
    cli::cli_abort(c(
      "{.arg s2_ind} must be 1-based: it indexes {.arg sd_est}, and R indexing starts at 1.",
      i = "Use {.code s2_ind = rep(1:2, each = n)} rather than {.code rep(0:1, each = n)}.",
      i = "Python impala uses 0-based indices; add 1 when porting a script."
    ))
  }
  if (min(s2_ind) < 1 || max(s2_ind) > ns2) {
    cli::cli_abort(c(
      "{.arg s2_ind} must lie between 1 and {.code length(sd_est)} = {ns2}.",
      i = "Observed range: {min(s2_ind)} to {max(s2_ind)}."
    ))
  }

  # Fill in the open interval when bounds are not supplied, so downstream
  # sampling never needs to branch on whether the user asked for bounds.
  if (is.null(sd_lower)) {
    sd_lower = rep(0, ns2)
  }
  if (is.null(sd_upper)) {
    sd_upper = rep(Inf, ns2)
  }
  if (length(sd_lower) != ns2 || length(sd_upper) != ns2) {
    cli::cli_abort(c(
      "{.arg sd_lower} and {.arg sd_upper} must have one entry per variance group.",
      i = "Expected length {ns2}, got {length(sd_lower)} and {length(sd_upper)}."
    ))
  }
  if (!is.numeric(sd_lower) || !is.numeric(sd_upper) ||
      anyNA(sd_lower) || anyNA(sd_upper)) {
    cli::cli_abort("{.arg sd_lower} and {.arg sd_upper} must be numeric with no missing values.")
  }
  if (any(sd_lower < 0)) {
    cli::cli_abort("{.arg sd_lower} must be non-negative.")
  }
  if (any(sd_lower >= sd_upper)) {
    cli::cli_abort("{.arg sd_lower} must be strictly less than {.arg sd_upper} in every group.")
  }
  # An out-of-bounds start leaves the M-H chain outside its own support, where
  # every candidate is rejected and the variance never moves.
  if (any(sd_est < sd_lower) || any(sd_est > sd_upper)) {
    cli::cli_abort(c(
      "{.arg sd_est} must lie within {.arg sd_lower} and {.arg sd_upper}.",
      i = "The sampler starts at {.arg sd_est}, so an out-of-bounds start cannot move."
    ))
  }

  if (is.null(theta_ind)) {
    theta_ind = rep(0, length(yobs))
  }

  vec = rep(0, ns2)
  for (i in 1:length(vec)) {
    vec[i] = sum(s2_ind == i)
  }
  if (any(vec == 0)) {
    empty = which(vec == 0)
    cli::cli_warn(c(
      "Variance group{?s} {empty} {?has/have} no observations in {.arg s2_ind}.",
      i = "Such groups are sampled from the prior alone."
    ))
  }

  if (is.null(D)) {
    nd = 0
  } else {
    nd = ncol(D)
  }

  model$exp_ind = theta_ind
  model$yobs = yobs
  if (!is.null(meas_error_cor)) {
    model$meas_error_cor = meas_error_cor
  }

  if (!is.null(D)) {
    model$D = D
    model$nd = nd
    model$discrep_tau = discrep_tau
  }


  obj$ntheta = c(obj$ntheta, length(unique(theta_ind)))
  obj$nclustmax = max(sum(obj$ntheta), 10)

  if (N == 0) {
    obj$ys = list(yobs)
    obj$y_lens = length(yobs)
    obj$theta_ind = list(theta_ind)
    obj$models = list(model)
    obj$sd_est = list(sd_est)
    obj$sd_lower = list(sd_lower)
    obj$sd_upper = list(sd_upper)
    obj$s2_df = list(s2_df)
    obj$ig_a = list(s2_df / 2)
    obj$ig_b = list(s2_df / 2 * sd_est^2)
    obj$s2_ind = list(s2_ind)
    obj$s2_exp_ind = list(1:ns2)
    obj$ns2 = list(ns2)
    obj$ny_s2 = list(vec)
    if (any(s2_df == 0)) {
      obj$s2_prior_kern = list(ldhc_kern)
    } else {
      obj$s2_prior_kern = list(ldig_kern)
    }
  } else {
    obj$ys[[N + 1]] = yobs
    obj$y_lens[[N + 1]] = length(yobs)
    obj$theta_ind[[N + 1]] = theta_ind
    obj$models[[N + 1]] = model
    obj$sd_est[[N + 1]] = sd_est
    obj$sd_lower[[N + 1]] = sd_lower
    obj$sd_upper[[N + 1]] = sd_upper
    obj$s2_df[[N + 1]] = s2_df
    obj$ig_a[[N + 1]] = s2_df / 2
    obj$ig_b[[N + 1]] = s2_df / 2 * sd_est^2
    obj$s2_ind[[N + 1]] = s2_ind
    obj$s2_exp_ind[[N + 1]] = 1:ns2
    obj$ns2[[N + 1]] = ns2
    obj$ny_s2[[N + 1]] = vec
    if (any(s2_df == 0)) {
      obj$s2_prior_kern[[N + 1]] = ldhc_kern
    } else {
      obj$s2_prior_kern[[N + 1]] = ldig_kern
    }
  }

  obj$nexp = obj$nexp + 1

  obj

}
