library(einsum)

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
  
  class(out) <- "AMcov_pool"
  out
}


update.AMcov_pool <- function(obj, x, m) {
  if (m > obj$start_adapt_iter) {
    obj$mu = obj$mu + (x[m - 1, , ] - obj$mu) / m
    tmp = x[m - 1, , ] - obj$mu
    obj$cov = ((m - 1) / m) * obj$cov  + ((m - 1) / (m * m)) * einsum('ti,tj->tij', tmp , tmp)
    obj$S   = obj$AM_SCALAR * einsum('ijk,i->ijk', obj$cov + diag(obj$p) * obj$eps, exp(obj$tau))
  } else if (m == obj$start_adapt_iter) {
    obj$mu = colMeans(x[1:m, , ])
    obj$cov = cov_3d_pcm(x[1:m, , ], obj$mu)
    obj$S   = obj$AM_SCALAR * einsum('ijk,i->ijk', obj$cov + diag(obj$p) * obj$eps, exp(obj$tau))
  }
  obj
}


update_tau.AMcov_pool <- function(obj, m) {
  if ((mod(m, 100) == 0) & (m > obj$start_adapt_iter)) {
    delta = min(0.5, 5 / sqrt(m + 1))
    obj$tau[obj$count_100 < 23] = obj$tau[obj$count_100 < 23] - delta
    obj$tau[obj$count_100 > 23] = obj$tau[obj$count_100 > 23] + delta
    obj$count_100 = obj$count_100 * 0
  }
  obj
}


gen_cand.AMcov_pool <- function(obj, x, m) {
  tmp = matrix(rnorm(obj$ntemps * obj$p), obj$ntemps)
  x_cand = x[m - 1, , ] + np.einsum('ijk,ik->ij', chol(obj$S), tmp)
  x_cand
}
