#' This function setups up temperature ladder for tempering
#' 
#' @param obj: CalibSetup Object
#' @param temperature_ladder: array of temperatures
#' @param start_temper: what MCMC sample to start tempering (default: 1000)
#' 
#' @return An object of class `CalibSetup`
#' 
setTemperatureLadder <- function(obj, temperature_ladder, start_temper) {
	obj$tl = temperature_ladder
	obj$itl = 1/obj$tl
	obj$ntemps = length(obj$tl)
	obj.nswap_per = floor(obj$ntemps/2)
	obj.start_temper = start_temper

	obj
}
