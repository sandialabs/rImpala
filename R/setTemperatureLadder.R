#' @title Set Temperature Ladder
#' @description Sets up the temperature ladder and swap scheme for parallel
#'   tempering.
#'
#' @param obj CalibSetup Object
#' @param temperature_ladder increasing vector of temperatures, the first equal
#'   to 1, e.g. `1.05^(0:49)`
#' @param start_temper what MCMC sample to start tempering (default: 1000)
#' @param swap_scheme how swap pairs are chosen in each swap round. `"deo"`
#'   (default) uses deterministic even/odd swaps between adjacent temperatures
#'   (non-reversible parallel tempering, Syed et al. 2022); `"seo"` picks the
#'   even/odd parity at random (the reversible counterpart); `"random"`
#'   reproduces the scheme of earlier versions (random disjoint pairs from the
#'   whole ladder) and gives the same chains for a fixed seed.
#' @param adapt_ladder if `TRUE`, re-space the temperatures during burn-in so
#'   that the swap rejection rate is the same between every pair of
#'   neighbouring temperatures (Syed et al. 2022, Algorithm 4). The coldest and
#'   hottest temperatures are kept fixed. Experimental (default: `FALSE`)
#' @param adapt_ladder_until last MCMC iteration at which the ladder may
#'   change. Defaults to `nburn` from [setMCMC()] if it is positive, otherwise
#'   `nmcmc \%/\% 2`. Discard samples up to this iteration as burn-in.
#'
#' @details With `swap_scheme = "deo"` and nearly every swap accepted (a very
#'   dense ladder), each replica returns to the cold chain every
#'   `2 * ntemps / nswap` iterations, which shows up as a spike in the cold
#'   chain's autocorrelation at that lag. It is not a bias, but it can alias
#'   with thinning; a sparser ladder removes it.
#'
#'   The samplers return the ladder actually used (`temperature_ladder`, after
#'   any adaptation), the estimated swap rejection rate on each adjacent pair
#'   (`swap_rejection`) and the number of cold-hot-cold round trips made by the
#'   replicas (`round_trips`).
#'
#' @references Syed, S., Bouchard-Cote, A., Deligiannidis, G. and Doucet, A.
#'   (2022). Non-reversible parallel tempering: a scalable highly parallel MCMC
#'   scheme. *Journal of the Royal Statistical Society Series B* 84(2),
#'   321-350.
#'
#'   Okabe, T., Kawata, M., Okamoto, Y. and Mikami, M. (2001). Replica-exchange
#'   Monte Carlo method for the isobaric-isothermal ensemble. *Chemical Physics
#'   Letters* 335, 435-439.
#'
#' @return An object of class `CalibSetup`
#'
#' @export
#'
setTemperatureLadder <- function(obj, temperature_ladder, start_temper = 1000,
                                 swap_scheme = c("deo", "seo", "random"),
                                 adapt_ladder = FALSE,
                                 adapt_ladder_until = NULL) {
  swap_scheme = match.arg(swap_scheme)
  if (isTRUE(adapt_ladder) && swap_scheme == "random") {
    cli::cli_abort("{.code adapt_ladder = TRUE} needs adjacent swaps; use {.code swap_scheme = \"deo\"} or {.val seo}.")
  }
  obj$tl = as.numeric(temperature_ladder)
  obj$itl = 1 / obj$tl
  obj$ntemps = length(obj$tl)
  obj$nswap_per = floor(obj$ntemps / 2)
  obj$start_temper = start_temper
  obj$swap_scheme = swap_scheme
  obj$adapt_ladder = isTRUE(adapt_ladder)
  obj$adapt_ladder_until = adapt_ladder_until

  obj
}
