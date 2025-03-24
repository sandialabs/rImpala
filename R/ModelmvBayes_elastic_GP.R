#' @title mvBayes elastic Emulator for Functional Outputs (can use different
#' BASS/BPPR type emulators)
#'
#' @description mvBayes_elastic Handles larger-dimensional functional responses (e.g., on
#' large spatial fields) using various inversion tricks. We require any
#' other covariance e.g., from discrepancy, measurement error, and basis
#' truncation error) to be diagonal. Contains a change on the likelihood
#' and built for a GP emulator. This function setups up emulator object.
#'
#' @param bmod a object of the type `mvBayes` of aligned functions
#' @param bmod_warp a object of the type `mvBayes` of warping functions
#' @param input_names cell array of strings of input variable names
#' @param exp_ind experiment indices (default: NULL)
#' @param s2 how to sample error variance (default: 'MH')
#' @param h h representation of warping function (default: FALSE)
#'
#' @return An object of class `ModelmvBayes_elastic_GP`
#'
#' @export
#'
ModelmvBayes_elastic_GP <- function(bmod,
                                    bmod_warp,
                                    input_names,
                                    exp_ind = NULL,
                                    s2 = 'MH',
                                    h = FALSE) {
  npc = bmod$basisInfo$nBasis
  is_mvBayes_available <- requireNamespace("mvBayes", quietly = TRUE)
  if (!is_mvBayes_available) {
    stop('install mvBayes for this basis option')
  }
  is_fdasrvf_available <- requireNamespace("fdasrvf", quietly = TRUE)
  if (!is_fdasrvf_available) {
    stop('install fdasrvf for this basis option')
  }

  if (s2 == 'gibbs') {
    cli::cli_abort("Cannot use Gibbs s2 for emulator models.")
  }

  if (is.null(exp_ind)) {
    exp_ind = 1
  }

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
    model_warp = bmod_warp,
    stochastic = TRUE,
    input_names = input_names,
    npc = npc,
    basis = bmod$basisInfo$basis,
    meas_error_corr = diag(nrow(bmod$basisInfo$basis)),
    discrep_cov = diag(nrow(bmod$basisInfo$basis)) * 1e-12,
    trunc_error_var = stats::cov(bmod$basisInfo$truncError) + stats::cov(bmod_warp$basisInfo$truncError),
    emu_vars = emu_vars,
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
    h = FALSE
  )

  class(obj) <- "ModelmvBayes_elastic_GP"

  obj
}


#' @export
step_m.ModelmvBayes_elastic_GP <- function(obj, ...) {
  obj
}


#' @export
discrep_sample.ModelmvBayes_elastic_GP <- function(obj, yobs, pred, cov, itemp, ...) {
  S = diag(obj$nd) / obj$discrep_tau + t(obj$D) %*% cov$inv %*% obj$D
  m = t(obj$D) %*% cov$inv %*% (yobs - pred)
  discrep_vars = chol_sample(solve(S, m), S / itemp)
  discrep_vars
}


#' @export
evalm.ModelmvBayes_elastic_GP <- function(obj,
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
    predf = stats::predict(obj$model, parmat_array)
    predv = stats::predict(obj$model_warp, parmat_array)

    if (dim(predf)[2] == 1) {
      if (obj$h) {
        gam = fdasrvf::h_to_gam(predv[1, , ])
      } else{
        gam = fdasrvf::v_to_gam(predv[1, , ])
      }

      M = dim(predf)[3]
      pred = fdasrvf::warp_f_gamma(predf[1, , ],
                                   seq(0, 1, length.out = M),
                                   fdasrvf::invertGamma(gam))

    } else {
      if (obj$h) {
        gam = fdasrvf::h_to_gam(t(predv[1, , ]))
      } else{
        gam = fdasrvf::v_to_gam(t(predv[1, , ]))
      }

      pred = predf[1, , ]

      for (i in 1:ncol(gam)) {
        pred[i, ] = fdasrvf::warp_f_gamma(predf[1, i, ],
                                          seq(0, 1, length.out = nrow(gam)),
                                          fdasrvf::invertGamma(gam[, i]))
      }
    }


  } else{
    cli::cli_abort("Not Implemented")
  }

  pred
}


#' @export
llik.ModelmvBayes_elastic_GP <- function(obj, yobs, pred, cov, ...) {
  vec = c(yobs - pred)
  out = -0.5 * (cov$ldet + t(vec) %*% cov$inv %*% vec)
  out
}


#' @export
lik_cov_inv.ModelmvBayes_elastic_GP <- function(obj, s2vec, ...) {
  N = length(s2vec)
  Sigma = cor2cov(obj$meas_error_corr, sqrt(s2vec))
  mat = Sigma + obj$trunc_error_var + obj$discrep_cov + obj$basis %*% diag(obj$emu_vars, nrow =
                                                                             obj$npc) %*% t(obj$basis)
  chol = chol(mat)
  ldet = 2 * sum(log(diag(chol)))
  inv = chol2inv(mat)
  out = list(inv = inv, ldet = ldet)
  out
}
