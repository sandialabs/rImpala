AMcov_pool <- function(ntemps,
                       p,
                       start_var,
                       start_adapt_iter,
                       tau_start) {
  S = array(0, dim = c(ntemps, p, p))
  for (i in 1:ntemps) {
    S[i, , ] = diag(p) * start_var
  }

  obj <- list(
    eps = 1e-12,
    AM_SCALAR = 2.4 ^ 2 / p,
    tau = rep(tau_start, each = ntemps),
    S = S,
    cov = array(0, dim = c(ntemps, p, p)),
    mu = matrix(0, ntemps, p),
    ntemps = ntemps,
    p = p,
    start_adapt_iter = start_adapt_iter,
    count_100 = rep(0, ntemps)
  )

  class(obj) <- "AMcov_pool"
  obj
}


#' @export
update_m.AMcov_pool <- function(obj, x, m, ...) {
  if (m > obj$start_adapt_iter) {
    obj$mu = obj$mu + (x[m - 1, , ] - obj$mu) / m
    tmp = x[m - 1, , ] - obj$mu
    if (ndims(tmp) == 0){
      tmp = matrix(tmp)
      a = ((m - 1) / (m * m)) * einsum::einsum('ti,tj->tij', tmp , tmp)
      obj$cov = ((m - 1) / m) * obj$cov  + matrix(a[,,1])
    } else {
      obj$cov = (((m-1)/m) * obj$cov  + ((m-1)/(m*m)) * einsum::einsum('ti,tj->tij', tmp , tmp))
    }

    eyetmp = replicate(dim(obj$cov)[1], diag(obj$p), simplify="array")
    if (obj$p == 1){
      eyetmp = matrix(eyetmp)
      tmp = array(obj$cov + eyetmp * obj$eps, dim=c(nrow(eyetmp),1,1))
      obj$S  = obj$AM_SCALAR * einsum::einsum('ijk,i->ijk', tmp, exp(obj$tau))
    } else {
      eyetmp = aperm(eyetmp, c(3,1,2))
      obj$S   = obj$AM_SCALAR * einsum::einsum('ijk,i->ijk', obj$cov + eyetmp * obj$eps, exp(obj$tau))
    }

  } else if (m == obj$start_adapt_iter) {
    obj$mu = colMeans(x[1:m, , ])
    obj$cov = cov_3d_pcm(x[1:m, , ], obj$mu)
    eyetmp = replicate(obj$ntemps, diag(obj$p), simplify="array")
    if (obj$p == 1){
      eyetmp = matrix(eyetmp)
      tmp = array(obj$cov + eyetmp * obj$eps, dim=c(length(eyetmp),1,1))
      obj$S  = obj$AM_SCALAR * einsum::einsum('ijk,i->ijk', tmp, exp(obj$tau))
    } else {
      eyetmp = aperm(eyetmp, c(3,1,2))
      obj$S   = obj$AM_SCALAR * einsum::einsum('ijk,i->ijk', obj$cov + eyetmp * obj$eps, exp(obj$tau))
    }
  }
  obj
}


#' @export
update_tau.AMcov_pool <- function(obj, m, ...) {
  if ((m %% 100 == 0) & (m > obj$start_adapt_iter)) {
    delta = min(0.5, 5 / sqrt(m + 1))
    obj$tau[obj$count_100 < 23] = obj$tau[obj$count_100 < 23] - delta
    obj$tau[obj$count_100 > 23] = obj$tau[obj$count_100 > 23] + delta
    obj$count_100 = obj$count_100 * 0
  }
  obj
}


#' @export
gen_cand.AMcov_pool <- function(obj, x, m, ...) {
  tmpchol = array(0, dim = dim(obj$S))
  for (i in 1:dim(obj$S)[1]) {
    tmpchol[i, , ] = t(chol(obj$S[i, , ]))
  }
  tmp = matrix(rnorm(obj$ntemps * obj$p), obj$ntemps)
  x_cand = x[m - 1, , ] + einsum::einsum('ijk,ik->ij', tmpchol, tmp)
  x_cand
}
