#' mvBayes Emulator for Functional Outputs (can use different
#' BASS/BPPR type emulators)
#'
#' ModelmvBayes Handles larger-dimensional functional responses (e.g., on
#' large spatial fields) using various inversion tricks. We require any
#' other covariance e.g., from discrepancy, measurement error, and basis
#' truncation error) to be diagonal
#'
#' This function setups up emulator object.
#'
#' @param bmod: a object of the type `mvBayes`
#' @param input_names: cell array of strings of input variable names
#' @param exp_ind: experiment indices (default: NULL)
#' @param s2: how to sample error variance (default: 'MH')
#' @param h: h representation of warping function (default: FALSE)
#'
#' @return An object of class `ModelmvBayes`
#'
#' @export
#'
ModelmvBayes <- function(bmod,
                         input_names,
                         exp_ind = NULL,
                         s2 = 'MH',
                         h = FALSE) {
  npc = bmod$basisInfo$nBasis
  if (class(bmod$bmList[[i]])=="bppr"){
    nmcmc = length(bmod$bmList[[1]]$sd_resid)
  } else {
    nmcmc = length(bmod$bmList[[1]]$s2)
  }


  if (s2 == 'gibbs') {
    cli::cli_abort("Cannot use Gibbs s2 for emulator models.")
  }

  if (is.null(exp_ind)) {
    exp_ind = 1
  }

  mod_s2 = matrix(0, nrow = nmcmc, npc)
  for (i in 1:npc) {
    if (class(bmod$bmList[[i]])=="bppr"){
      mod_s2[, i] = bmod$bmList[[i]]$sd_resid^2
    } else {
      mod_s2[, i] = bmod$bmList[[i]]$s2
    }

  }

  obj <- list(
    model = bmod,
    stochastic = TRUE,
    nmcmc = nmcmc,
    input_names = input_names,
    basis = bmod$basisInfo$basis,
    meas_error_corr = diag(nrow(bmod$basisInfo$basis)),
    discrep_cov = diag(nrow(bmod$basisInfo$basis)) * 1e-12,
    ii = 1,
    trunc_error_var = cov(bmod$basisInfo$truncError),
    mod_s2 = mod_s2,
    emu_vars = mod_s2[1, ],
    yobs = NULL,
    marg_lik_cov = NULL,
    discrep_vars = NULL,
    nd = 0,
    discrep_tau = 1,
    D = NULL,
    discrep = 0,
    exp_ind = exp_ind,
    nexp = max(exp_ind),
    s2 = s2
  )

  class(obj) <- "ModelmvBayes"

  obj
}


#' @export
step_m.ModelmvBayes <- function(obj) {
  obj$ii = sample(obj$nmcmc, 1)
  obj$emu_vars = obj$mod_s2[obj$ii, ]
  obj
}


#' @export
discrep_sample.ModelmvBayes <- function(obj, yobs, pred, cov, itemp) {
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% cov$inv %*% D
  m = t(obj$D) %*% cov$inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S,m), S / itemp)
  discrep_vars
}


#' @export
evalm.ModelmvBayes <- function(obj,
                               parmat,
                               pool = TRUE,
                               nugget = FALSE) {
  fn = obj$input_names
  parmat_array = matrix(0, length(parmat[[fn[1]]]), length(fn))
  for (i in 1:length(fn)) {
    parmat_array[, i] = parmat[[fn[i]]]
  }

  if (pool) {
    pred = predict(obj$model,
                   parmat_array,
                   mcmc.use = obj$ii,
                   nugget = nugget)
  } else{
    cli::cli_abort("Not Implemented")
  }

  pred[1, , ]
}


#' @export
llik.ModelmvBayes <- function(obj, yobs, pred, cov) {
  vec = c(yobs - pred)
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelmvBayes <- function(obj, s2vec) {
  N = length(s2vec)
  Sigma = cor2cov(obj$meas_error_corr, sqrt(s2vec))
  mat = Sigma + obj$trunc_error_var + obj$discrep_cov + obj$basis %*% diag(obj$emu_vars) %*% t(obj$basis)
  chol = chol(mat)
  ldet = 2 * sum(log(diag(chol)))
  inv = chol2inv(chol)
  out = list(inv = inv, ldet = ldet)
  out
}
