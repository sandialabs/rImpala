#' PCA Based Model Emulator using BASS
#'
#' This function setups up emulator object
#'
#' @param bmod: a object of the type `bassBasis`
#' @param input_names: cell array of strings of input variable names
#' @param exp_ind: experiment indices (default: NULL)
#' @param s2: how to sample error variance (default: 'MH')
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
    trunc_error_var = diag(cov(t(bmod$dat$trunc.error)))
  } else {
    # need to check this
    trunc_error_var = diag(cov(t(bmod$dat$trunc.error)))[1]
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
    s2 = s2
  )

  class(obj) <- "ModelBassPca_func"

  obj
}


#' @export
step_m.ModelBassPca_func <- function(obj) {
  obj$ii = sample(obj$nmcmc, 1)
  obj$emu_vars = obj$mod_s2[obj$ii, ]
  obj
}


#' @export
discrep_sample.ModelBassPca_func <- function(obj, yobs, pred, cov, itemp) {
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% cov$inv %*% D
  m = t(obj$D) %*% cov$inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S) %*% m, S / itemp)
  discrep_vars
}


#' @export
evalm.ModelBassPca_func <- function(obj,
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
llik.ModelBassPca_func <- function(obj, yobs, pred, cov) {
  vec = c(yobs - pred)
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelBassPca_func <- function(obj, s2vec) {
  vec = obj$trunc_error_var + s2vec
  Ainv = diag(1 / vec)
  Aldet = sum(log(vec))
  out = swm(Ainv,
            obj$basis,
            diag(1 / obj$emu_vars),
            t(obj$basis),
            Aldet,
            sum(log(obj$emu_vars)))
  out
}
