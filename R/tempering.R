# Swap scheduling and temperature-ladder adaptation for parallel tempering,
# shared by calibPool(), calibHier() and calibClust(). Port of
# src/impala/superCal/tempering.py from python impala.
#
# Swap schemes
#   "deo" (default): deterministic even/odd swaps between adjacent
#     temperatures, the non-reversible parallel tempering of Okabe et al.
#     (2001) analysed by Syed et al. (2022). Rounds alternate between the pairs
#     (1,2), (3,4), ... and (2,3), (4,5), .... Because the parity alternates
#     deterministically, an accepted state keeps moving in the same direction
#     along the ladder, so the ladder is crossed in O(N) swap rounds instead of
#     the O(N^2) of a reversible random walk.
#   "seo": stochastic even/odd, the reversible counterpart of "deo".
#   "random": the scheme of earlier versions: each round draws
#     floor(ntemps / 2) disjoint pairs uniformly from the whole ladder.
#
# Ladder adaptation
#   The rejection rate r_k between inverse temperatures beta_k and beta_{k+1}
#   is estimated with the Rao-Blackwellised estimator 1 - min(1, alpha_k). The
#   cumulative rejection Lambda(beta) (the communication barrier) is then made
#   equally spaced along a new ladder (Syed et al. 2022, Algorithm 4). The
#   coldest and hottest inverse temperatures stay fixed. Adaptation happens at
#   the end of rounds whose length doubles, starting at 100 iterations, and
#   stops at `adapt_until`, after which the sampler is a fixed Markov kernel.
#
# References
#   Okabe, T., Kawata, M., Okamoto, Y. and Mikami, M. (2001). Chemical Physics
#     Letters 335, 435-439.
#   Syed, S., Bouchard-Cote, A., Deligiannidis, G. and Doucet, A. (2022).
#     Non-reversible parallel tempering: a scalable highly parallel MCMC scheme.
#     JRSS-B 84(2), 321-350.

SWAP_SCHEMES <- c("deo", "seo", "random")


# Disjoint adjacent pairs (k, k + 1), 1-based, with (k - 1) %% 2 == parity, as
# a 2 x npair matrix (row 1 the colder temperature).
even_odd_pairs <- function(ntemps, parity) {
  if (parity + 1 > ntemps - 1) {
    return(matrix(integer(0), 2, 0))
  }
  first = seq.int(parity + 1, ntemps - 1, by = 2)
  rbind(first, first + 1L, deparse.level = 0)
}


# Re-space a decreasing inverse-temperature ladder (itl[1] = 1) so that the
# estimated cumulative rejection is equally spaced along it. `rejection` has
# length(itl) - 1 entries, one per adjacent edge. End points are kept.
communication_barrier_ladder <- function(itl, rejection, min_rejection = 1e-4) {
  itl = as.numeric(itl)
  rejection = pmin(pmax(as.numeric(rejection), min_rejection), 1)
  # increasing grid beta_0 < ... < beta_N (hot -> cold)
  beta = rev(itl)
  rej = rev(rejection)
  barrier = c(0, cumsum(rej))
  targets = seq(0, barrier[length(barrier)], length.out = length(beta))
  new_beta = stats::approx(barrier, beta, xout = targets, ties = "ordered")$y
  new_beta[1] = beta[1]
  new_beta[length(new_beta)] = beta[length(beta)]
  rev(new_beta)
}


# Create the swap scheduler. It is an environment, so the samplers update it in
# place. The samplers must use `sched$itl` (a working copy of the ladder) rather
# than `setup$itl`, and rebuild anything derived from it whenever
# sched_maybe_adapt() returns TRUE.
#
# `random_layout` reproduces the pair layout each sampler used before swap
# schemes existed, so swap_scheme = "random" gives the same chain as before for
# a fixed seed: "split" pairs x[s] with x[s + npair] (calibPool), "interleaved"
# pairs x[2s - 1] with x[2s] (calibHier, calibClust).
make_swap_scheduler <- function(setup, random_layout = c("interleaved", "split"),
                                first_round = 100) {
  random_layout = match.arg(random_layout)
  scheme = setup$swap_scheme
  if (is.null(scheme)) scheme = "deo"
  adapt = isTRUE(setup$adapt_ladder)
  if (!scheme %in% SWAP_SCHEMES) {
    cli::cli_abort("{.arg swap_scheme} must be one of {.val {SWAP_SCHEMES}}, not {.val {scheme}}.")
  }
  if (adapt && scheme == "random") {
    cli::cli_abort("Ladder adaptation needs adjacent swaps; use {.code swap_scheme = \"deo\"} or {.val seo}.")
  }
  until = setup$adapt_ladder_until
  if (is.null(until)) {
    nburn = if (is.null(setup$nburn)) 0 else setup$nburn
    until = if (nburn > 0) nburn else setup$nmcmc %/% 2
  }
  ntemps = length(setup$itl)
  nedge = max(ntemps - 1, 0)

  s = new.env(parent = emptyenv())
  s$itl = as.numeric(setup$itl)
  s$ntemps = ntemps
  s$scheme = scheme
  s$random_layout = random_layout
  s$adapt = adapt && ntemps > 2
  s$adapt_until = until
  s$round_length = first_round
  s$next_adapt = NULL
  s$n_rounds = 0
  s$n_adaptations = 0
  s$rej_sum = rep(0, nedge)
  s$rej_n = rep(0, nedge)
  # statistics accumulated since the last adaptation (or since the start)
  s$rejection_total = rep(0, nedge)
  s$attempts_total = rep(0, nedge)
  # round trips: which replica sits at each temperature, and whether it has
  # reached the hottest temperature since it was last at the coldest
  s$label = seq_len(ntemps)
  s$seen_hot = rep(FALSE, ntemps)
  if (ntemps > 1) s$seen_hot[ntemps] = TRUE
  s$round_trips = 0
  s
}


# Pairs of temperature indices for the next swap round: a 2 x npair matrix,
# one pair per column, row 1 the colder temperature for adjacent schemes.
sched_pairs <- function(s) {
  s$n_rounds = s$n_rounds + 1
  if (s$scheme == "random") {
    npair = s$ntemps %/% 2
    x = sample(1:s$ntemps, 2 * npair, FALSE)
    if (s$random_layout == "split") {
      return(t(matrix(x, npair, 2)))
    }
    return(matrix(x, 2))
  }
  parity = if (s$scheme == "deo") (s$n_rounds - 1) %% 2 else sample.int(2, 1) - 1
  even_odd_pairs(s$ntemps, parity)
}


# Record one round of swap attempts: `log_alpha` the log acceptance ratio of
# each pair (column of `pairs`) and `accepted` which of them were accepted.
sched_record <- function(s, pairs, log_alpha, accepted) {
  if (ncol(pairs) == 0) {
    return(invisible(s))
  }
  if (any(accepted)) {
    sched_track_round_trips(s, pairs[, accepted, drop = FALSE])
  }
  lo = pmin(pairs[1, ], pairs[2, ])
  adjacent = abs(pairs[2, ] - pairs[1, ]) == 1
  if (!any(adjacent)) {
    return(invisible(s))
  }
  edge = lo[adjacent]
  la = as.numeric(log_alpha)[adjacent]
  la[is.na(la)] = -Inf
  rej = 1 - exp(pmin(la, 0))
  for (k in seq_along(edge)) {
    e = edge[k]
    s$rej_sum[e] = s$rej_sum[e] + rej[k]
    s$rej_n[e] = s$rej_n[e] + 1
    s$rejection_total[e] = s$rejection_total[e] + rej[k]
    s$attempts_total[e] = s$attempts_total[e] + 1
  }
  invisible(s)
}


sched_track_round_trips <- function(s, swapped) {
  for (k in seq_len(ncol(swapped))) {
    a = swapped[1, k]
    b = swapped[2, k]
    s$label[c(a, b)] = s$label[c(b, a)]
  }
  s$seen_hot[s$label[s$ntemps]] = TRUE
  cold = s$label[1]
  if (s$seen_hot[cold]) {
    s$round_trips = s$round_trips + 1
    s$seen_hot[cold] = FALSE
  }
  invisible(s)
}


# Estimated rejection rate on every adjacent edge (NaN if never attempted).
sched_rejection_rates <- function(s) {
  s$rejection_total / s$attempts_total
}


# Adapt the ladder if an adaptation round has just finished at iteration m.
# Returns TRUE if s$itl changed.
sched_maybe_adapt <- function(s, m) {
  if (!s$adapt) {
    return(FALSE)
  }
  if (!is.null(s$adapt_until) && m > s$adapt_until) {
    return(FALSE)
  }
  if (is.null(s$next_adapt)) {
    s$next_adapt = m + s$round_length
    return(FALSE)
  }
  if (m < s$next_adapt || any(s$rej_n == 0)) {
    return(FALSE)
  }
  rej = s$rej_sum / s$rej_n
  s$itl = communication_barrier_ladder(s$itl, rej)
  s$rej_sum[] = 0
  s$rej_n[] = 0
  # report post-adaptation statistics only
  s$rejection_total[] = 0
  s$attempts_total[] = 0
  s$round_length = s$round_length * 2
  s$next_adapt = m + s$round_length
  s$n_adaptations = s$n_adaptations + 1
  TRUE
}


# Log prior kernel of s2 for experiment i to use in tempering swap ratios. It
# must be the prior the within-chain s2 update actually targets. The conjugate
# Gibbs update draws from the inverse-gamma full conditional with shape ig_a and
# scale ig_b (both 0 when s2_df == 0, i.e. p(s2) proportional to 1/s2), so the
# swap uses the inverse-gamma kernel with those parameters. The Metropolis
# update uses setup$s2_prior_kern[[i]] directly (the half-Cauchy-type kernel
# when any s2_df == 0), so the swap uses that too.
s2_swap_kern <- function(setup, i) {
  if (identical(setup$models[[i]]$s2, "gibbs")) {
    return(ldig_kern)
  }
  setup$s2_prior_kern[[i]]
}
