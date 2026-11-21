#' @title Set Hierarchical Priors
#'
#' @description Sets the hyperpriors used by [calibHier()] and [calibClust()].
#'   Each parameter vector \eqn{\theta_j} (one per value of `theta_ind`, see
#'   [addVecExperiments()]) is drawn from
#'   \eqn{N(\theta_0, \Sigma_0)}, with
#'   \eqn{\theta_0 \sim N(}`theta0_prior_mean`, `theta0_prior_cov`\eqn{)} and
#'   \eqn{\Sigma_0 \sim} Inverse-Wishart(`Sigma0_prior_df`,
#'   `Sigma0_prior_scale`). All of these live on the 0-1 scale the sampler
#'   works on, not the native scale of `bounds`.
#'
#' @param obj `CalibSetup` object
#' @param theta0_prior_mean prior mean of \eqn{\theta_0}, length `p`
#'   (default: `rep(0.5, p)`)
#' @param theta0_prior_cov prior covariance of \eqn{\theta_0}, `p x p`
#'   (default: `diag(p) * 100`, effectively flat on the unit cube)
#' @param Sigma0_prior_df Inverse-Wishart degrees of freedom for
#'   \eqn{\Sigma_0}; must exceed `p + 1` so the prior mean exists. Larger values
#'   pull the \eqn{\theta_j} together (default: `p + 20`)
#' @param Sigma0_prior_scale Inverse-Wishart scale matrix, `p x p`. The prior
#'   mean of \eqn{\Sigma_0} is `Sigma0_prior_scale / (Sigma0_prior_df - p - 1)`
#'   (default: `diag(p) * 0.1^2`)
#' @param theta0_start optional `ntemps x p` matrix of 0-1 scale starting values
#'   for \eqn{\theta_0}, one row per temperature (default: `NULL`, a random
#'   start satisfying the constraints)
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
setHierPriors <- function(obj,
                          theta0_prior_mean = NULL,
                          theta0_prior_cov = NULL,
                          Sigma0_prior_df = NULL,
                          Sigma0_prior_scale = NULL,
                          theta0_start = NULL) {
  p = obj$p
  if (is.null(theta0_prior_mean)) theta0_prior_mean = rep(0.5, p)
  if (is.null(theta0_prior_cov)) theta0_prior_cov = diag(p) * 100
  if (is.null(Sigma0_prior_df)) Sigma0_prior_df = p + 20
  if (is.null(Sigma0_prior_scale)) Sigma0_prior_scale = diag(p) * 0.1^2

  theta0_prior_cov = as.matrix(theta0_prior_cov)
  Sigma0_prior_scale = as.matrix(Sigma0_prior_scale)

  if (length(theta0_prior_mean) != p) {
    cli::cli_abort("{.arg theta0_prior_mean} must have length {p}, one entry per parameter.")
  }
  if (!all(dim(theta0_prior_cov) == p) || !all(dim(Sigma0_prior_scale) == p)) {
    cli::cli_abort("{.arg theta0_prior_cov} and {.arg Sigma0_prior_scale} must be {p} x {p} matrices.")
  }
  if (!is_pos_def(theta0_prior_cov) || !is_pos_def(Sigma0_prior_scale)) {
    cli::cli_abort("{.arg theta0_prior_cov} and {.arg Sigma0_prior_scale} must be symmetric positive definite.")
  }
  if (length(Sigma0_prior_df) != 1 || !is.finite(Sigma0_prior_df) ||
      Sigma0_prior_df <= p + 1) {
    cli::cli_abort(c(
      "{.arg Sigma0_prior_df} must be a single number greater than {.code p + 1} = {p + 1}.",
      i = "The sampler starts Sigma0 at its prior mean, which needs {.code Sigma0_prior_df > p + 1}."
    ))
  }
  if (!is.null(theta0_start)) {
    theta0_start = as.matrix(theta0_start)
    if (ncol(theta0_start) == 1 && p > 1) theta0_start = t(theta0_start)
    if (ncol(theta0_start) != p || any(theta0_start < 0) || any(theta0_start > 1)) {
      cli::cli_abort(c(
        "{.arg theta0_start} must be an {.code ntemps x p} matrix on the 0-1 scale.",
        i = "A single row is recycled across temperatures."
      ))
    }
  }

  obj$theta0_prior_mean = as.numeric(theta0_prior_mean)
  obj$theta0_prior_cov = theta0_prior_cov
  obj$Sigma0_prior_df = Sigma0_prior_df
  obj$Sigma0_prior_scale = Sigma0_prior_scale
  obj$theta0_start = theta0_start

  obj
}


#' @title Set Cluster Priors
#'
#' @description Sets the Dirichlet process settings used by [calibClust()].
#'   The \eqn{\theta_j} are grouped into at most `nclustmax` clusters that share
#'   a parameter value. The concentration parameter \eqn{\eta} has a
#'   Gamma(`eta_prior_shape`, `eta_prior_rate`) prior; larger \eqn{\eta} favours
#'   more clusters. Call [setHierPriors()] as well, since the cluster values are
#'   drawn from \eqn{N(\theta_0, \Sigma_0)}.
#'
#' @param obj `CalibSetup` object
#' @param nclustmax maximum number of clusters (default: `NULL`, meaning
#'   `max(sum(ntheta), 10)`, where `ntheta` counts the parameter vectors set up
#'   through `theta_ind`)
#' @param eta_prior_shape shape of the Gamma prior on \eqn{\eta} (default: 2)
#' @param eta_prior_rate rate of the Gamma prior on \eqn{\eta} (default: 0.1)
#' @param scale_proposal if `TRUE`, divide each cluster's proposal covariance by
#'   its number of members, which can help when cluster sizes are very uneven
#'   (default: `FALSE`)
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
setClusterPriors <- function(obj,
                             nclustmax = NULL,
                             eta_prior_shape = 2,
                             eta_prior_rate = 0.1,
                             scale_proposal = FALSE) {
  if (is.null(nclustmax)) {
    nclustmax = max(sum(obj$ntheta), 10)
  }
  if (length(nclustmax) != 1 || nclustmax < 1 || nclustmax != as.integer(nclustmax)) {
    cli::cli_abort("{.arg nclustmax} must be a positive whole number.")
  }
  if (eta_prior_shape <= 0 || eta_prior_rate <= 0) {
    cli::cli_abort("{.arg eta_prior_shape} and {.arg eta_prior_rate} must be positive.")
  }
  obj$nclustmax = as.integer(nclustmax)
  obj$eta_prior_shape = eta_prior_shape
  obj$eta_prior_rate = eta_prior_rate
  obj$scale_proposal = isTRUE(scale_proposal)
  obj
}


is_pos_def <- function(x) {
  isSymmetric(unname(x), tol = 1e-8) &&
    !inherits(try(chol(x), silent = TRUE), "try-error")
}
