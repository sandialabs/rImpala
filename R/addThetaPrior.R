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
  if(is.null(pname) | !(pname %in% pnames)){
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
#' @examples
#' # Joint normal prior on t_1 and t_2
#' my_joint_prior <- function(params) {
#'   # params is a list like list(t_1=0.5, t_2=0.3)
#'   mvtnorm::dmvnorm(c(params$t_1, params$t_2),
#'                    mean=c(0, 0),
#'                    sigma=matrix(c(1, 0.5, 0.5, 1), 2, 2),
#'                    log=TRUE)
#' }
#' setup <- addJointThetaPrior(setup, c("t_1", "t_2"), my_joint_prior)
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
