#' mvBayes Emulator for Functional Outputs (GP)
#'
#' ModelmvBayes_GP Handles larger-dimensional functional responses (e.g., on
#' large spatial fields) using various inversion tricks. We require any
#' other covariance e.g., from discrepancy, measurement error, and basis
#' truncation error) to be diagonal. Built for a GP emulator.
#'
#' This function setups up emulator object.
#'
#' @param bmod a object of the type `mvBayes`
#' @param input_names cell array of strings of input variable names
#' @param exp_ind experiment indices (default: NULL)
#' @param s2 how to sample error variance (default: 'MH')
#' @param h h representation of warping function (default: FALSE)
#'
#' @return An object of class `ModelmvBayes_GP`
#'
#' @export
#'
ModelmvBayes_GP <- function(bmod,
                            input_names,
                            exp_ind = NULL,
                            s2 = 'MH',
                            h = FALSE) {
  if (s2 == 'gibbs') {
    cli::cli_abort("Cannot use Gibbs s2 for emulator models.")
  }
  is_mvBayes_available <- requireNamespace("mvBayes", quietly = TRUE)
  if (!is_mvBayes_available) {
    stop('install mvBayes for this basis option')
  }

  if (is.null(exp_ind)) {
    exp_ind = 1
  }

  npc = bmod$basisInfo$nBasis
  emu_vars = rep(NA, npc)
  N = length(bmod$bmList[[1]]$covparms)
  for (ii in 1:npc) {
    tmp = stats::predict(bmod$bmList[[ii]],
                         bmod$bmList[[ii]]$locs,
                         joint = FALSE,
                         predvar = TRUE)
    emu_vars[ii] = mean(tmp$vars)
  }

  obj <- list(
    model = bmod,
    stochastic = TRUE,
    input_names = input_names,
    npc = npc,
    basis = t(bmod$basisInfo$basis),
    meas_error_corr = diag(ncol(bmod$basisInfo$basis)),
    discrep_cov = diag(ncol(bmod$basisInfo$basis)) * 1e-12,
    trunc_error_var = stats::cov(bmod$basisInfo$truncError),
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

  class(obj) <- "ModelmvBayes_GP"

  obj
}


#' @export
step_m.ModelmvBayes_GP <- function(obj, ...) {
  obj
}


#' @export
discrep_sample.ModelmvBayes_GP <- function(obj, yobs, pred, cov, itemp, ...) {
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% cov$inv %*% obj$D
  m = t(obj$D) %*% cov$inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S, m), S / itemp)
  discrep_vars
}


#' @export
evalm.ModelmvBayes_GP <- function(obj,
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
    pred = stats::predict(obj$model, parmat_array)
  } else{
    cli::cli_abort("Not Implemented")
  }

  if (ndims(pred) == 3){
    return(pred[1, , ])
  } else {
    return(pred)
  }
}


#' @export
llik.ModelmvBayes_GP <- function(obj, yobs, pred, cov, ...) {
  vec = c(yobs - pred)
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelmvBayes_GP <- function(obj, s2vec, ...) {
  N = length(s2vec)
  Sigma = cor2cov(obj$meas_error_corr, sqrt(s2vec))
  mat = Sigma + obj$trunc_error_var + obj$discrep_cov + obj$basis %*% diag(obj$emu_vars, nrow =
                                                                             obj$npc) %*% t(obj$basis)
  chol = chol(mat)
  ldet = 2 * sum(log(diag(chol)))
  inv = chol2inv(chol)
  out = list(inv = inv, ldet = ldet)
  out
}
