library(einsum)

ndims <- function(x){
  return(length(dim(x)))
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
  mtmp = bounds[, 1]
  diff = (bounds[, 2] - bounds[, 1])
  out = (x - mtmp) / diff
  out
}


unnormalize <- function(z, bounds) {
  mtmp = bounds[, 1]
  diff = (bounds[, 2] - bounds[, 1])
  out = z * diff + mtmp
  out
}


tran_unif <- function(th, bounds, names) {
  out = list()
  tbounds = unnormalize(th, bounds)
  for (i in 1:length(names)) {
    out[[names[i]]] = tbounds[, i]
  }
  out
}


cov_3d_pcm <- function(arr, mean) {
  N = nrow(arr)
  if (ndims(arr) == 3){
    meantmp = replicate(N, mean, simplify="array")
    meantmp = aperm(meantmp, c(3,1,2))
    out = einsum('kij,kil->ijl', arr - meantmp, arr - meantmp) / (N - 1)
  } else if (ndims(arr) == 2){
    meantmp = replicate(N, mean, simplify="array")
    meantmp = aperm(meantmp, c(2,1))
    out = matrix(diag(cov(arr - meantmp)))
  }
  out
}
