#' @title Set MCMC Parameters
#' @description This function setups up MCMC parameters for adaptive MCMC,
#' also includes tempering and decorrelation steps
#'
#' @param obj `calibPool` object
#' @param nmcmc number of mcmc iterations
#' @param nburn number of mcmc burn in iterations (default: 0)
#' @param thin number of samples to thin (default: 1)
#' @param decor number of mcmc iterations before decorrelation step (default: 100)
#' @param start_var_theta start variance of theta proposal (default: 1e-8)
#' @param start_tau_theta start tau of theta (default: 0)
#' @param start_var_ls2 start variance of sigma proposal (default: 1e-5)
#' @param start_tau_ls2 start tau of sigma (default: 0)
#' @param start_adapt_iter number of iterations before adaption
#'                          (default: 300)
#' @return An object of class `CalibSetup`
#'
#' @export
#'
setMCMC <- function(obj,
                    nmcmc,
                    nburn = 0,
                    thin = 1,
                    decor = 100,
                    start_var_theta = 1e-8,
                    start_tau_theta = 0,
                    start_var_ls2 = 1e-5,
                    start_tau_ls2 = 0,
                    start_adapt_iter = 300) {
  obj$nmcmc = nmcmc
  obj$nburn = nburn
  obj$thin = thin
  obj$decor = decor
  obj$start_var_theta = start_var_theta
  obj$start_tau_theta = start_tau_theta
  obj$start_var_ls2 = start_var_ls2
  obj$start_tau_ls2 = start_tau_ls2
  obj$start_adapt_iter = start_adapt_iter

  obj
}
