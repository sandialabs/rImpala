llik <- function(obj, ...) {
  UseMethod("llik")
}


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


step <- function(obj, ...) {
  UseMethod("step")
}


step.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


eval <- function(obj, ...) {
  UseMethod("eval")
}


eval.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


update <- function(obj, ...){
	UseMethod("update")
}


update.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


update_tau <- function(obj, ...){
	UseMethod("update_tau")
}


update_tau.default <- function(obj, ...) {
  cat("This is a generic function\n")
}


gen_cand <- function(obj, ...){
	UseMethod("gen_cand")
}


gen_cand.default <- function(obj, ...) {
  cat("This is a generic function\n")
}
