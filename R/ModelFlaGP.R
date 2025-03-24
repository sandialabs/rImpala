#' @title FlaGP Emulator for Functional Outputs (GP)
#'
#' @description ModelFlaGP Handles larger-dimensional functional responses (e.g., on
#' large spatial fields) using various inversion tricks. We require any
#' other covariance e.g., from discrepancy, measurement error, and basis
#' truncation error) to be diagonal. This function setups up emulator object.
#'
#' @param bmod a object of the type `flagp`
#' @param input_names cell array of strings of input variable names
#' @param exp_ind experiment indices (default: NULL)
#' @param s2 how to sample error variance (default: 'MH')
#' @param h h representation of warping function (default: FALSE)
#'
#' @return An object of class `ModelFlaGP`
#'
#' @export
#'
ModelFlaGP <- function(bmod,
                       input_names,
                       exp_ind = NULL,
                       s2 = 'MH',
                       h = FALSE) {

  if (s2 == 'gibbs') {
    cli::cli_abort("Cannot use Gibbs s2 for emulator models.")
  }

  if (is.null(exp_ind)) {
    exp_ind = 1
  }

  npc = bmod$basis$sim$n.pc
  X.orig = bmod$XT.data$sim$X$orig
  tmp1 = FLaGP::predict(bmod,X.pred.orig=X.orig,verbose=F)
  emu_vars = stats::var(t(tmp1$y.mean - bmod$Y.data$sim$orig))
  truncError = bmod$basis$sim$B %*% bmod$basis$sim$V.t -  bmod$Y.data$sim$trans
  trunc_error_var = stats::cov(t(truncError))

  obj <- list(
    model = bmod,
    stochastic = TRUE,
    input_names = input_names,
    basis = t(bmod$basis$sim$B),
    meas_error_corr = diag(nrow(bmod$basis$sim$B)),
    discrep_cov = diag(nrow(bmod$basis$sim$B)) * 1e-12,
    trunc_error_var = trunc_error_var,
    yobs = NULL,
    marg_lik_cov = NULL,
    discrep_vars = NULL,
    nd = 0,
    ii = 1,
    emu_vars = emu_vars,
    discrep_tau = 1,
    D = NULL,
    discrep = 0,
    exp_ind = exp_ind,
    nexp = max(exp_ind),
    s2 = s2
  )

  class(obj) <- "ModelFlaGP"

  obj
}


#' @export
step_m.ModelFlaGP <- function(obj, ...) {
  obj
}


#' @export
discrep_sample.ModelFlaGP <- function(obj, yobs, pred, cov, itemp, ...) {
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% cov$inv %*% obj$D
  m = t(obj$D) %*% cov$inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S,m), S / itemp)
  discrep_vars
}


#' @export
evalm.ModelFlaGP <- function(obj,
                             parmat,
                             pool = TRUE,
                             nugget = FALSE,
                             ...) {
  fn = obj$input_names
  parmat_array = matrix(0, length(parmat[[fn[1]]]), length(fn))
  for (i in 1:length(fn)) {
    parmat_array[, i] = parmat[[fn[i]]]
  }

  if (pool) {
    pred = predict(obj$model,
                   X.pred.orig=parmat_array,
                   verbose=F)$y.mean
  } else{
    cli::cli_abort("Not Implemented")
  }

  t(pred)
}


#' @export
llik.ModelFlaGP <- function(obj, yobs, pred, cov, ...) {
  vec = c(yobs - pred)
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelFlaGP <- function(obj, s2vec, ...) {
  N = length(s2vec)
  Sigma = cor2cov(obj$meas_error_corr, sqrt(s2vec))
  mat = Sigma + obj$trunc_error_var + obj$discrep_cov + obj$emu_vars
  chol = chol(mat)
  ldet = 2 * sum(log(diag(chol)))
  inv = chol2inv(chol)
  out = list(inv = inv, ldet = ldet)
  out
}
