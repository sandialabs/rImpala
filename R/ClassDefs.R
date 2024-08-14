llik <- function(obj, ...) {
  UseMethod("llik")
}


#' @export
llik.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


lik_cov_inv <- function(obj, ...) {
  UseMethod("lik_cov_inv")
}


lik_cov_inv.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


discrep_sample <- function(obj, ...) {
  UseMethod("discrep_sample")
}


discrep_sample.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


step_m <- function(obj, ...) {
  UseMethod("step_m")
}

#' @export
step_m.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


#' @export
evalm <- function(obj, ...) {
  UseMethod("evalm")
}


#' @export
evalm.default <- function(obj, ...) {
  cat("This is a generic function\n")
}

update_m <- function(obj, ...) {
  UseMethod("update_m")
}

#' @export
update_m.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


update_tau <- function(obj, ...) {
  UseMethod("update_tau")
}


update_tau.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


gen_cand <- function(obj, ...) {
  UseMethod("gen_cand")
}


gen_cand.default <- function(obj, ...) {
  cat("This is a generic function\n")
}
