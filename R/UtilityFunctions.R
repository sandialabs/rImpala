ndims <- function(x) {
  return(length(dim(x)))
}


cor2cov <- function(V, sd) {
  V * tcrossprod(sd)
}


ldhc_kern <- function(x, a, b) {
  out = -log(x + 1)
  out
}


ldig_kern <- function(x, a, b) {
  out = (-a - 1) * log(x) - b / x
  out
}


swm <- function(Ainv, U, Cinv, V, Aldet, Cldet) {
  in_mat = chol_solve(Cinv + V %*% Ainv %*% U)
  inv1 = Ainv - Ainv %*% U %*% in_mat$inv %*% V %*% Ainv
  ldet = in_mat$ldet + Aldet + Cldet
  out = list(inv = inv1, ldet = ldet)
  out
}


chol_solve <- function(x) {
  R = chol(x)
  ldet = 2 * sum(log(diag(R)))
  inv1 = solve(x)
  out = list(inv = inv1, ldet = ldet)
  out
}


chol_sample <- function(mean, cov) {
  out = mean + t(chol(cov)) %*% stats::rnorm(length(mean))
  return(out)
}


#' @title Compare to bounds
#' @description This function compares variables to bounds
#'
#' @param x list of of parameters
#' @param bounds list of bounds for each parameter that is a two parameter vector with high and low
#'
#' @return a vector of `TRUE` or `FALSE` if the values are within the bounds
#'
#' @export
#'
cf_bounds <- function(x, bounds) {
  k = names(bounds)
  good = x[[k[1]]] < bounds[[k[1]]][2]
  for (i in 1:length(k)) {
    good = good * (x[[k[i]]] < bounds[[k[i]]][2]) * (x[[k[i]]] > bounds[[k[i]]][1])
  }
  good = as.logical(good)

  good
}


normalize <- function(x, bounds) {
  m = bounds[, 1]
  diff = (bounds[, 2] - bounds[, 1])
  diff_tmp = t(replicate(nrow(x), diff))
  mtmp = t(replicate(nrow(x), m))
  out = (x - mtmp) / diff
  out
}


unnormalize <- function(z, bounds) {
  m = bounds[, 1]
  diff = (bounds[, 2] - bounds[, 1])
  diff_tmp = t(replicate(nrow(z), diff))
  if (nrow(diff_tmp) == 1 & nrow(z) > 1){
    diff_tmp = t(diff_tmp)
  }
  mtmp = t(replicate(nrow(z), m))
  if (nrow(mtmp) == 1 & nrow(z) > 1){
    mtmp = t(mtmp)
  }
  out = z * diff_tmp + mtmp
  out
}


#' @title Transform parameters
#' @description Transforms variables from the scale 0 to 1 to variable provided by bounds
#'
#' @param th list of of parameters
#' @param bounds list of bounds for each parameter that is a two parameter vector with high and low
#' @param names vector of variable names
#'
#' @return a list of parameters
#'
#' @export
#'
tran_unif <- function(th, bounds, names) {
  out = list()
  tbounds = unnormalize(th, bounds)
  for (i in 1:length(names)) {
    out[[names[i]]] = tbounds[, i]
  }
  out
}


# Evaluate the log-prior for calibration parameters. `theta` is a named list of
# parameter vectors (one entry per temperature), as returned by `tran_unif`.
# `priors` is `setup$theta_prior`; when NULL/empty this returns a vector of zeros
# so calibration is unaffected when no prior has been added.
eval_theta_priors = function(theta, priors, tnames=NULL){
  lp = rep(0,length(theta[[1]]))
  if(!is.null(tnames))
    names(theta) = tnames

  for(p in priors){

    # Check if this is a joint prior (has 'names' field) or independent prior (has 'name' field)
    if(!is.null(p$names)){
      # Joint prior - extract multiple parameters and call custom function
      params_list = lapply(p$names, function(nm) theta[[nm]])
      names(params_list) = p$names
      add = p$log_density_fn(params_list)
    } else {
      # Independent prior - use standard distributions
      x = theta[[p$name]]

      add <- switch(
        p$dist,
        normal   = stats::dnorm(x, mean = p$params$mean, sd = p$params$sd, log = TRUE),
        lognormal  = stats::dlnorm(x, meanlog = p$params$meanlog, sdlog = p$params$sdlog, log = TRUE),
        beta   = stats::dbeta(x, shape1 = p$params$shape1, shape2 = p$params$shape2, log = TRUE),
        uniform   = stats::dunif(x, min = p$params$min, max = p$params$max, log = TRUE),
        gamma  = stats::dgamma(x, shape = p$params$shape, rate = p$params$rate, log = TRUE),
        cauchy = stats::dcauchy(x, location = p$params$location, scale = p$params$scale, log = TRUE),
        stop("Unsupported dist: ", p$dist)
      )
    }

    lp = lp + add
  }
  lp
}


cov_3d_pcm <- function(arr, mean) {
  N = nrow(arr)
  if (ndims(arr) == 3) {
    meantmp = replicate(N, mean, simplify = "array")
    meantmp = aperm(meantmp, c(3, 1, 2))
    out = einsum::einsum('kij,kil->ijl', arr - meantmp, arr - meantmp) / (N - 1)
  } else if (ndims(arr) == 2) {
    meantmp = replicate(N, mean, simplify = "array")
    meantmp = aperm(meantmp, c(2, 1))
    tmp = array(arr - meantmp, dim = c(nrow(arr), ncol(arr), 1))
    out = einsum::einsum('kij,kil->ijl', tmp, tmp) / (N - 1)
    out = out[, , 1]
  }
  out
}
