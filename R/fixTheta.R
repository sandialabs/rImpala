#' @title Fix calibration parameters
#'
#' @description This method holds one or more calibration parameters at fixed
#'   values, so that [calibPool()] samples the conditional posterior of the
#'   remaining parameters. This supports conditional inference and cut-Bayes
#'   style analyses, where some inputs are treated as known, or as determined by
#'   a separate module, rather than learned from the data at hand.
#'
#' @param obj `CalibSetup` Object
#' @param pname character vector of parameter names to fix; each must appear in
#'               `names(obj$bounds)`
#' @param value numeric vector of values to fix the parameters at, on the native
#'               scale of `obj$bounds` and the same length as `pname`
#'
#' @return An object of class `CalibSetup`
#'
#' @details Fixed parameters keep their column in the `theta` and `theta_native`
#'   output of [calibPool()]; that column is simply constant at the fixed value.
#'   The adaptive Metropolis proposal and the decorrelation step both act on the
#'   free parameters only, so fixing a parameter does not degrade the mixing of
#'   the others.
#'
#'   Calling `fixTheta()` again for a parameter that is already fixed replaces
#'   its value. Use [unfixTheta()] to return a parameter to being sampled.
#'
#'   A prior added with [addThetaPrior()] or [addJointThetaPrior()] for a
#'   parameter that is fixed contributes a constant to the log posterior and so
#'   has no effect on the samples; a warning is issued when that is detected.
#'
#'   Fixing every parameter is allowed: [calibPool()] then samples the error
#'   variances and any discrepancy terms with `theta` held at the fixed values.
#'
#' @export
#'
#' @examples
#' bounds = list()
#' bounds[['t_1']] = c(0, 1)
#' bounds[['t_2']] = c(-2, 2)
#' setup <- CalibSetup(bounds, cf_bounds)
#'
#' # hold t_2 at 0.5 and calibrate t_1 conditional on it
#' setup <- fixTheta(setup, "t_2", 0.5)
#' setup$theta_fixed
#'
#' # several at once
#' setup <- fixTheta(setup, c("t_1", "t_2"), c(0.25, -1))
#' setup$theta_fixed
#'
#' # and back to sampling both
#' setup <- unfixTheta(setup)
#' setup$theta_fixed
#'
fixTheta <- function(obj, pname, value) {
  pnames = names(obj$bounds)

  if (missing(pname) ||
      !is.character(pname) || length(pname) == 0) {
    stop('pname must be a character vector of parameter names')
  }
  if (!all(pname %in% pnames)) {
    stop(
      'Parameter not in set of input names for any model in setup$models: ',
      paste(setdiff(pname, pnames), collapse = ', ')
    )
  }
  if (anyDuplicated(pname) > 0) {
    stop('pname contains duplicate parameter names')
  }
  if (missing(value) ||
      !is.numeric(value) || length(value) != length(pname)) {
    stop('value must be a numeric vector the same length as pname')
  }
  if (any(!is.finite(value))) {
    stop('value must be finite')
  }

  idx = match(pname, pnames)
  outside = (value < obj$bounds_mat[idx, 1]) |
    (value > obj$bounds_mat[idx, 2])
  if (any(outside)) {
    stop('Fixed value outside of bounds for: ',
         paste(pname[outside], collapse = ', '))
  }

  fixed = obj$theta_fixed
  if (is.null(fixed)) {
    fixed = stats::setNames(numeric(0), character(0))
  }
  fixed[pname] = value
  # keep the fixed set in bounds order so downstream indexing is predictable
  obj$theta_fixed = fixed[pnames[pnames %in% names(fixed)]]

  warn_fixed_prior(obj, pname)

  if (length(obj$theta_fixed) == obj$p) {
    cli::cli_warn(
      paste0(
        "All ",
        obj$p,
        " calibration parameters are fixed. calibPool() will sample the ",
        "error variances and any discrepancy terms only."
      )
    )
  }

  obj
}


#' @title Unfix calibration parameters
#'
#' @description This method releases parameters fixed by [fixTheta()] so that
#'   [calibPool()] samples them again.
#'
#' @param obj `CalibSetup` Object
#' @param pname character vector of parameter names to release, or `NULL`
#'               (the default) to release every fixed parameter
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
unfixTheta <- function(obj, pname = NULL) {
  if (is.null(pname)) {
    obj$theta_fixed = NULL
    return(obj)
  }

  if (!is.character(pname) || length(pname) == 0) {
    stop('pname must be a character vector of parameter names')
  }
  if (!all(pname %in% names(obj$bounds))) {
    stop(
      'Parameter not in set of input names for any model in setup$models: ',
      paste(setdiff(pname, names(obj$bounds)), collapse = ', ')
    )
  }

  not_fixed = setdiff(pname, names(obj$theta_fixed))
  if (length(not_fixed) > 0) {
    cli::cli_warn(paste0(
      "Parameter(s) not currently fixed: ",
      paste(not_fixed, collapse = ', ')
    ))
  }

  keep = setdiff(names(obj$theta_fixed), pname)
  if (length(keep) == 0) {
    obj$theta_fixed = NULL
  } else {
    obj$theta_fixed = obj$theta_fixed[keep]
  }

  obj
}


# Warn when a prior has been attached to a parameter that is now fixed: the
# prior then contributes a constant to the log posterior and cancels out of
# every acceptance ratio, so it silently does nothing. `pname` limits the check
# to the parameters just fixed.
warn_fixed_prior <- function(obj, pname) {
  if (is.null(obj$theta_prior)) {
    return(invisible(NULL))
  }

  prior_names = unlist(lapply(obj$theta_prior, function(p) {
    if (is.null(p$names)) p$name else p$names
  }))
  hit = intersect(pname, prior_names)

  if (length(hit) > 0) {
    cli::cli_warn(
      paste0(
        "Prior set for fixed parameter(s): ",
        paste(hit, collapse = ', '),
        ". A prior on a fixed parameter is constant and does not affect the ",
        "posterior."
      )
    )
  }

  invisible(NULL)
}


# Resolve the fixed/free split of the calibration parameters for a CalibSetup.
# Returns the positions of the fixed and free parameters within
# `names(setup$bounds)` plus the fixed values mapped onto the 0-1 scale
# `calibPool` samples on. Setups built before `theta_fixed` existed have no such
# field, so a NULL there means nothing is fixed.
fixed_theta_split <- function(setup) {
  pnames = names(setup$bounds)
  fixed = setup$theta_fixed
  unit = rep(NA_real_, setup$p)
  names(unit) = pnames

  if (length(fixed) > 0) {
    idx = match(names(fixed), pnames)
    lo = setup$bounds_mat[idx, 1]
    hi = setup$bounds_mat[idx, 2]
    unit[idx] = (fixed - lo) / (hi - lo)
  }

  # unname the indices: they are used to subscript matrix columns, where stray
  # names would propagate into the output
  list(fixed_idx = unname(which(!is.na(unit))),
       free_idx = unname(which(is.na(unit))),
       unit = unit[!is.na(unit)])
}


# Overwrite the fixed columns of an (n x p) matrix of unit-scale parameters with
# their fixed values. A no-op when nothing is fixed.
set_fixed_cols <- function(x, split) {
  if (length(split$fixed_idx) > 0) {
    x[, split$fixed_idx] = matrix(split$unit,
                                  nrow(x),
                                  length(split$fixed_idx),
                                  byrow = TRUE)
  }
  x
}
