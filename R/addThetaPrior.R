#' @title Add theta prior
#'
#' @description This method adds a prior for the specified theta
#'
#' @param obj `CalibSetup` Object
#' @param dist string defining distribution e.g., 'normal', 'beta', etc.
#' @param params list of parameters e.g., if dist='normal', params is a list with elements 'mean' and 'sd'
#' @param pname name of parameter
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
addThetaPrior <- function(obj,
                          dist='uniform',
                          params=list(min=0,max=1),
                          pname=NULL
                          ) {
  pnames = names(obj$bounds)
  if(is.null(pname) || length(pname) != 1 || !(pname %in% pnames)){
    stop('No parameter name given or parameter not in set of input names for any model in setup$models')
  }

  if(is.null(obj$theta_prior))
    nprior = 0
  else
    nprior = length(obj$theta_prior)

  if (nprior == 0) {
    obj$theta_prior = list(list(name=pname,dist=dist,params=params))
  } else {
    obj$theta_prior[[nprior + 1]] = list(name=pname,dist=dist,params=params)
  }

  obj

}


#' @title Add joint theta prior
#'
#' @description This method adds a joint prior over multiple theta parameters
#'
#' @param obj `CalibSetup` Object
#' @param pnames character vector of parameter names to apply joint prior to
#' @param log_density_fn function that takes named list of parameter values and returns log-density
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
#' @details
#' `log_density_fn` receives a named list of parameter *vectors*, one element per
#' temperature, and must return a vector of log densities of that same length.
#' Write it vectorised, or wrap a scalar density in [mapply()].
#'
#' @examples
#' # Joint prior on t_1 and t_2 correlating them through a bivariate normal.
#' # params arrives as list(t_1 = <vector>, t_2 = <vector>), one entry per
#' # temperature, so the result must be a vector of the same length.
#' my_joint_prior <- function(params) {
#'   rho <- 0.5
#'   z <- (params$t_1^2 - 2 * rho * params$t_1 * params$t_2 + params$t_2^2) /
#'     (1 - rho^2)
#'   -0.5 * z - log(2 * pi * sqrt(1 - rho^2))
#' }
#'
#' bounds = list()
#' bounds[['t_1']] = c(0, 1)
#' bounds[['t_2']] = c(0, 1)
#' setup <- CalibSetup(bounds, cf_bounds)
#' setup <- addJointThetaPrior(setup, c("t_1", "t_2"), my_joint_prior)
#'
#' # returns one log density per temperature
#' my_joint_prior(list(t_1 = c(0.2, 0.5), t_2 = c(0.3, 0.6)))
#'
addJointThetaPrior <- function(obj,
                               pnames,
                               log_density_fn) {

  # Validate parameter names
  valid_names = names(obj$bounds)
  if(!all(pnames %in% valid_names)){
    stop('One or more parameter names not in set of input names for models in setup')
  }

  # Validate function
  if(!is.function(log_density_fn)){
    stop('log_density_fn must be a function')
  }

  if(is.null(obj$theta_prior))
    nprior = 0
  else
    nprior = length(obj$theta_prior)

  if (nprior == 0) {
    obj$theta_prior = list(list(names=pnames, log_density_fn=log_density_fn))
  } else {
    obj$theta_prior[[nprior + 1]] = list(names=pnames, log_density_fn=log_density_fn)
  }

  obj
}
