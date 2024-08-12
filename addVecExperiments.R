#' This method adds vector experiments to calibration object
#'
#' @param obj: CalibSetup Object
#' @param yobs: a vector of the experiment or observation
#' @param model: emulator (currently expecting a object of class
#'               ModelBassPca_func, ModelBassPca_func_mf, or
#'               ModelBpprPca_func
#' @param sd_est: estimate of standard deviation
#' @param s2_df: degrees of freedom of inverse gamma prior
#' @param s2_ind: indices of function
#' @param meas_error_cor: measurement error correlation (default: NULL)
#' @param theta_ind: indices of theta (default: NULL)
#' @param D: discrepency basis (matrix of columns of basis, default:
#'           NaN)
#' @param discrep_tau: discrepency sampling tau
#'
#' @return An object of class `CalibSetup`
#'
addVecExperiments <- function(obj,
                              yobs,
                              model,
                              sd_est,
                              s2_df,
                              s2_ind,
                              meas_error_cor = NULL,
                              theta_ind = NULL,
                              D = NULL,
                              discrep_tau = 1) {
  N = length(obj$ys)
  
  if (is.null(theta_ind)) {
    theta_ind = rep(0, length(yobs))
  }
  
  vec = rep(0, length(sd_est))
  for (i in 1:length(vec)) {
    vec[i] = sum(s2_ind == i)
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
    obj$s2_df = list(s2_df)
    obj$ig_a = list(s2_df / 2)
    obj$ig_b = list(s2_df / 2 * sd_est ^ 2)
    obj$s2_ind = list(s2_ind)
    obj$s2_exp_ind = list(1:length(sd_est))
    obj$ns2 = list(length(sd_est))
    obj$ny_s2 = list(vec)
    if (sum(s2_df == 0) > 1) {
      obj.s2_prior_kern = list(ldhc_kern)
    } else {
      obj.s2_prior_kern = list(ldig_kern)
    }
  } else {
    obj$ys[[N + 1]] = yobs
    obj$y_lens[[N + 1]] = length(yobs)
    obj$theta_ind[[N + 1]] = theta_ind
    obj$models[[N + 1]] = model
    obj$sd_est[[N + 1]] = sd_est
    obj$s2_df[[N + 1]] = s2_df
    obj$ig_a[[N + 1]] = s2_df / 2
    obj$ig_b[[N + 1]] = s2_df / 2 * sd_est ^ 2
    obj$s2_ind[[N + 1]] = s2_ind
    obj$s2_exp_ind[[N + 1]] = 1:length(sd_est)
    obj$ns2[[N + 1]] = length(sd_est)
    obj$ny_s2[[N + 1]] = vec
    if (sum(s2_df == 0) > 1) {
      obj$s2_prior_kern[[N + 1]] = ldhc_kern
    } else {
      obj$s2_prior_kern[[N + 1]] = ldig_kern
    }
  }
  
  obj$nexp = obj$nexp + 1
  
  obj
  
}
