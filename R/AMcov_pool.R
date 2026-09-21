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
    AM_SCALAR = 2.4^2 / p,
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


# Rebuild the (ntemps x p x p) proposal covariance from the running covariance,
# ridged by `eps` and scaled per temperature by exp(tau).
scale_S <- function(obj) {
  eyetmp = array(0, dim = c(obj$ntemps, obj$p, obj$p))
  for (i in 1:obj$ntemps) {
    eyetmp[i, , ] = diag(obj$p)
  }
  obj$AM_SCALAR * einsum::einsum('ijk,i->ijk',
                                 obj$cov + eyetmp * obj$eps,
                                 exp(obj$tau))
}


#' @export
update_m.AMcov_pool <- function(obj, x, m, cols = NULL, ...) {
  # `cols` restricts the adaptation to a subset of the columns of `x`, used when
  # some calibration parameters are held fixed: the object tracks only the free
  # ones, so obj$p is length(cols) rather than dim(x)[3].
  if (is.null(cols)) {
    cols = seq_len(obj$p)
  }

  if (m > obj$start_adapt_iter) {
    # keep the (ntemps x p) shape: single-index slices of a 3-d array drop dims
    xprev = matrix(x[m - 1, , cols], obj$ntemps, obj$p)
    obj$mu = obj$mu + (xprev - obj$mu) / m
    tmp = xprev - obj$mu
    obj$cov = ((m - 1) / m) * obj$cov +
      ((m - 1) / (m * m)) * einsum::einsum('ti,tj->tij', tmp, tmp)
    obj$S = scale_S(obj)

  } else if (m == obj$start_adapt_iter) {
    xhist = array(x[1:m, , cols], dim = c(m, obj$ntemps, obj$p))
    obj$mu = matrix(colMeans(xhist), obj$ntemps, obj$p)
    obj$cov = array(cov_3d_pcm(xhist, obj$mu),
                    dim = c(obj$ntemps, obj$p, obj$p))
    obj$S = scale_S(obj)
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
gen_cand.AMcov_pool <- function(obj, x, m, cols = NULL, ...) {
  # as in update_m, `cols` selects the free columns of `x` when some calibration
  # parameters are fixed; the returned candidate has obj$p columns
  if (is.null(cols)) {
    cols = seq_len(obj$p)
  }

  tmpchol = array(0, dim = dim(obj$S))
  for (i in 1:dim(obj$S)[1]) {
    tmpchol[i, , ] = t(chol(obj$S[i, , ]))
  }
  tmp = matrix(stats::rnorm(obj$ntemps * obj$p), obj$ntemps)
  x_cand = matrix(x[m - 1, , cols], obj$ntemps, obj$p) +
    einsum::einsum('ijk,ik->ij', tmpchol, tmp)
  x_cand
}
