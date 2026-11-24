#' @title Pooled Bayesian Model Calibration
#'
#' @description This function setups up calibration object
#' @param bounds a list with fields of variable names with values of
#'                two dimensional arrays where first element is lower bound
#'                and second element is upper bound.
#' @param constraint_func function handle to constraint function on
#'                         variables, called as `f(x, bounds)` or, when it
#'                         takes a third argument, `f(x, bounds, constants)`
#'                         with the constants of the first model added (see
#'                         [constraints_ptw()]). Defaults to [cf_bounds()],
#'                         which only checks the bounds.
#'
#' @return An object of class `CalibSetup` which is a list with the following
#'   components:
#'
#' - `nexp`: number of experiments
#'
#' - `theta_fixed`: named vector of parameters held fixed during calibration,
#'   `NULL` until [fixTheta()] is called
#'
#' @export
#'

CalibSetup <- function(bounds, constraint_func = cf_bounds) {
  fn = names(bounds)
  bmattmp = matrix(0, length(bounds), 2)
  for (i in 1:length(bounds)) {
    bmattmp[i, ] = bounds[[fn[i]]]
  }

  out <- list(
    nexp = 0,
    ys = NULL,
    y_lens = NULL,
    models = NULL,
    tl = 1,
    itl = 1 / 1,
    bounds = bounds,
    bounds_mat = bmattmp,
    checkConstraints = constraint_func,
    constants = NULL,
    p = length(bounds),
    nmcmc = 10000,
    nburn = 5000,
    thin = 5,
    decor = 100,
    ntemps = 1,
    sd_est = NULL,
    sd_lower = NULL,
    sd_upper = NULL,
    s2_df = NULL,
    ig_a = NULL,
    ig_b = NULL,
    s2_ind = NULL,
    s2_exp_ind = NULL,
    ns2 = NULL,
    ny_s2 = NULL,
    ntheta = c(),
    theta_ind = NULL,
    nswap = 5,
    nswap_per = 0,
    start_temper = Inf,
    swap_scheme = "deo",
    adapt_ladder = FALSE,
    adapt_ladder_until = NULL,
    theta_prior = NULL,
    theta_fixed = NULL,
    s2_prior_kern = NULL,
    wt = NULL,
    # hierarchical / clustered calibration, see setHierPriors()/setClusterPriors()
    theta0_prior_mean = NULL,
    theta0_prior_cov = NULL,
    Sigma0_prior_df = NULL,
    Sigma0_prior_scale = NULL,
    theta0_start = NULL,
    nclustmax = NULL,
    eta_prior_shape = 2,
    eta_prior_rate = 0.1,
    scale_proposal = FALSE
  )

  class(out) <- "CalibSetup"

  out
}
