#' Pooled Bayesian Model Calibration
#' 
#' This function setups up calibration object
#' 
#' @param bounds: a list with fields of variable names with values of 
#'                two dimensional arrays where first element is lower bound 
#'                and second element is upper bound.
#' @param constraint_func: function handle to constraint function on 
#'                         variables. Defaults to function on checking 
#'                         bounds
#' 
#' @return An object of class `CalibSetup` which is a list with the following
#'   components:
#'
#' - `nexp`: number of experiments
#' 
#' @export

CalibSetup <- function(bounds, constraint_func){
	out <- list(
		nexp = 0,
		ys = NULL,
		y_lens = NULL,
		models = NULL,
		tl = 1,
		itl = 1/1,
		bounds = bounds,
		p = length(bounds),
		nmcmc = 10000,
		nburn = 5000,
		thin = 5,
		decor = 100,
		ntemps = 1,
		sd_est = NULL,
		s2_df = NULL,
		ig_a = NULL,
		ig_b = NULL,
		s2_ind = NULL,
		s2_exp_ind = NULL,
		ns2 = NULL,
		s2_exp_ind = NULL,
		ns2 = NULL,
		ny_s2 = NULL,
		ntheta = c(),
		theta_ind = NULL,
		nswap = 5,
		s2_prior_kern = NULL
  	)

	class(out) <- "CalibSetup"

  	out
}
