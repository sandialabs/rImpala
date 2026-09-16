#' @title Log likelihood
#' @description Generic for the log likelihood of an emulator model. Methods are
#'   provided for each supported emulator class; supply your own to calibrate
#'   with a custom emulator.
#'
#' @param obj an emulator model object
#' @param ... additional arguments passed to the method, typically the
#'   observations, the emulator prediction, and the marginal covariance returned
#'   by [lik_cov_inv()]
#'
#' @return the log likelihood, a scalar
#'
#' @export
#'
llik <- function(obj, ...) {
  UseMethod("llik")
}


#' @export
llik.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Likelihood covariance inverse
#' @description Generic that builds the inverse and log determinant of the
#'   marginal likelihood covariance for an emulator model.
#'
#' @param obj an emulator model object
#' @param ... additional arguments passed to the method, typically the vector of
#'   error variances
#'
#' @return a list with elements `inv` (the inverse covariance) and `ldet` (its
#'   log determinant)
#'
#' @export
#'
lik_cov_inv <- function(obj, ...) {
  UseMethod("lik_cov_inv")
}


#' @export
lik_cov_inv.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Sample discrepancy coefficients
#' @description Generic that draws discrepancy basis coefficients for an
#'   emulator model with a discrepancy basis attached.
#'
#' @param obj an emulator model object
#' @param ... additional arguments passed to the method, typically the
#'   observations, the prediction, the marginal covariance, and the inverse
#'   temperature
#'
#' @return a vector of sampled discrepancy coefficients
#'
#' @export
#'
discrep_sample <- function(obj, ...) {
  UseMethod("discrep_sample")
}


#' @export
discrep_sample.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Advance a stochastic emulator
#' @description Generic that advances an emulator one MCMC step, used for
#'   emulators that carry their own posterior samples.
#'
#' @param obj an emulator model object
#' @param ... additional arguments passed to the method
#'
#' @return the updated emulator model object
#'
#' @export
#'
step_m <- function(obj, ...) {
  UseMethod("step_m")
}

#' @export
step_m.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title evalm constructor
#' @description Default constructor for evalm class
#'
#' @param obj evalm object
#' @param ... additional arguments passed to method
#'
#' @return An object of class `evalm`
#'
#' @export
#'
evalm <- function(obj, ...) {
  UseMethod("evalm")
}


#' @export
evalm.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Update adaptive covariance
#' @description Generic that updates the running mean and covariance of an
#'   adaptive Metropolis proposal.
#'
#' @param obj an adaptive covariance object, e.g. from `AMcov_pool`
#' @param ... additional arguments passed to the method, typically the chain of
#'   samples and the current iteration
#'
#' @return the updated adaptive covariance object
#'
#' @export
#'
update_m <- function(obj, ...) {
  UseMethod("update_m")
}

#' @export
update_m.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Update proposal scaling
#' @description Generic that applies diminishing adaptation to the proposal
#'   scale based on the recent acceptance rate.
#'
#' @param obj an adaptive covariance object, e.g. from `AMcov_pool`
#' @param ... additional arguments passed to the method, typically the current
#'   iteration
#'
#' @return the updated adaptive covariance object
#'
#' @export
#'
update_tau <- function(obj, ...) {
  UseMethod("update_tau")
}


#' @export
update_tau.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}


#' @title Generate a proposal
#' @description Generic that draws a Metropolis candidate from an adaptive
#'   proposal covariance.
#'
#' @param obj an adaptive covariance object, e.g. from `AMcov_pool`
#' @param ... additional arguments passed to the method, typically the chain of
#'   samples and the current iteration
#'
#' @return a matrix of candidate values, one row per temperature
#'
#' @export
#'
gen_cand <- function(obj, ...) {
  UseMethod("gen_cand")
}


#' @export
gen_cand.default <- function(obj, ...) {
  cli::cli_alert_info("This is a generic function\n")
}
