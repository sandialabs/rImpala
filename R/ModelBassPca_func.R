#' @title PCA Based Model Emulator using BASS
#'
#' @description This function setups up emulator object
#' @param bmod a object of the type `bassBasis`
#' @param input_names cell array of strings of input variable names
#' @param exp_ind experiment indices (default: NULL)
#' @param s2 how to sample error variance (default: 'MH')
#'
#' @return An object of class `ModelBassPca_func`
#'
#' @export
#'
ModelBassPca_func <- function(bmod,
                              input_names,
                              exp_ind = NULL,
                              s2 = 'MH') {
  npc = ncol(bmod$dat$basis)
  nmcmc = length(bmod$mod.list[[1]]$s2)

  if (s2 == 'gibbs') {
    cli::cli_abort("Cannot use Gibbs s2 for emulator models.")
  }

  if (is.null(exp_ind)) {
    exp_ind = 1
  }

  if (npc > 1) {
    trunc_error_var = diag(stats::cov(t(bmod$dat$trunc.error)))
  } else {
    # need to check this
    trunc_error_var = diag(stats::cov(t(bmod$dat$trunc.error)))[1]
  }

  mod_s2 = matrix(0, nrow = nmcmc, npc)
  for (i in 1:npc) {
    mod_s2[, i] = bmod$mod.list[[i]]$s2
  }

  obj <- list(
    model = bmod,
    stochastic = TRUE,
    nmcmc = length(bmod$mod.list[[1]]$s2),
    input_names = input_names,
    npc = npc,
    basis = bmod$dat$basis,
    meas_error_corr = diag(nrow(bmod$dat$basis)),
    discrep_cov = diag(nrow(bmod$dat$basis)) * 1e-12,
    ii = 1,
    trunc_error_var = trunc_error_var,
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
    s2 = s2,
    supports_wt = TRUE
  )

  class(obj) <- "ModelBassPca_func"

  obj
}


#' @export
step_m.ModelBassPca_func <- function(obj, ...) {
  obj$ii = sample(obj$nmcmc, 1)
  obj$emu_vars = obj$mod_s2[obj$ii, ]
  obj
}


#' @export
discrep_sample.ModelBassPca_func <- function(obj, yobs, pred, cov, itemp, wt = NULL, ...) {
  inv = weight_inv(cov$inv, wt)
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% inv %*% obj$D
  m = t(obj$D) %*% inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S) %*% m, S / itemp)
  discrep_vars
}


#' @export
evalm.ModelBassPca_func <- function(obj,
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
    pred = stats::predict(obj$model,
                            parmat_array,
                            mcmc.use = obj$ii,
                            nugget = nugget)
  } else{
    cli::cli_abort("Not Implemented")
  }

  # keep the (ntemps x ny) shape: pred[1, , ] drops to a bare vector whenever
  # ntemps == 1 or ny == 1, which then fails downstream in calibPool
  matrix(pred[1, , ], nrow = dim(pred)[2], ncol = dim(pred)[3])
}


#' @export
llik.ModelBassPca_func <- function(obj, yobs, pred, cov, wt = NULL, ...) {
  # weights enter as a power likelihood: residuals scaled by sqrt(wt), and the
  # diagonal part of the log determinant by wt (see lik_cov_inv)
  vec = c(yobs - pred)
  if (!is.null(wt)) {
    vec = sqrt(wt) * vec
  }
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelBassPca_func <- function(obj, s2vec, wt = NULL, inds = NULL, ...) {
  # `inds` selects the observations s2vec belongs to, so hierarchical and
  # clustered calibration can build the covariance of one theta's subset
  if (is.null(inds)) {
    inds = seq_along(s2vec)
  }
  if (is.null(wt)) {
    wt = rep(1, length(s2vec))
  }
  tev = obj$trunc_error_var
  if (length(tev) > 1) {
    tev = tev[inds]
  }
  vec = tev + s2vec
  Ainv = diag(1 / vec, nrow = length(vec))
  Aldet = sum(wt * log(vec))
  basis = obj$basis[inds, , drop = FALSE]
  out = swm(
    Ainv,
    basis,
    diag(1 / obj$emu_vars, nrow = obj$npc),
    t(basis),
    Aldet,
    sum(log(obj$emu_vars))
  )
  out
}
