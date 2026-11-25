# impala (development version)

Ports lanl/impala#69 and #71 from python impala.

* parallel tempering now uses non-reversible deterministic even/odd swaps
  between adjacent temperatures by default (`setTemperatureLadder(swap_scheme =
  "deo")`); `"seo"` and `"random"` (the previous behaviour, same chains for a
  fixed seed) are also available. Optional ladder adaptation during burn-in
  (`adapt_ladder = TRUE`, experimental). `calibPool()`, `calibHier()` and
  `calibClust()` return `temperature_ladder`, `swap_rejection` and
  `round_trips`
* fix: the tempering swap ratio used the half-Cauchy s2 kernel when
  `s2_df = 0` even though the Gibbs s2 update targets p(s2) ∝ 1/s2, which
  biased the cold chain once swaps started (all samplers)
* fix: the Metropolis s2 update tempered the Jacobian of its random walk on
  log(s2), so hot chains targeted a different family than the swap ratio
  assumed and the cold chain was biased (about -6% in a test case; all
  samplers)
* fix: `calibClust()` now tempers only the likelihood (and s2 prior). It
  previously tempered the normal density of occupied clusters and the theta0 /
  Sigma0 updates but not empty clusters, memberships or the hyperpriors, so the
  chains had no common target and swaps biased the cold chain
* port the material strength physics models from python impala:
  `MaterialModel()` with PTW, Johnson-Cook, Steinberg-Guinan and constant flow
  stress plus specific heat, density, melt temperature and shear modulus
  sub-models (`getoptions_ModelMaterialStrength()` lists them), and
  `ModelMaterialStrength()` to calibrate them against stress-strain curves.
  `constraints_ptw()` checks the PTW parameter ordering
* add `ModelF()` to calibrate an R function directly, with a diagonal Gaussian
  likelihood
* add `hierTruncationWeights()`, which diagnoses the contraction that
  `calibHier()`'s truncation of the hierarchy to the bounds induces on
  `theta0` and `Sigma0`. The truncation is applied without its normalising
  constant, so the effective hyperprior is the stated one times
  `prod_j P(theta_j in B | theta0, Sigma0)`, which rewards a central `theta0`
  and a small `Sigma0`. The function reweights the saved cold-chain draws by
  that factor to the power `-sum(ntheta)` and reports `ess_frac`, `tilt_nats`
  and raw-versus-reweighted summaries, so a finished run can be corrected
  retrospectively where the reweighting is effective. `weightedQuantile()` is
  exported alongside it (ports lanl/impala#70)
* `CalibSetup()` now defaults `constraint_func` to `cf_bounds()`. A constraint
  function that takes a third argument receives the first model's constants,
  so it can check parameters held constant
* port the remaining python impala tests: physics model regressions and FLAG
  comparisons, pooled/hierarchical/clustered PTW calibration, the Friedman
  uniform-posterior check and the theta prior checks

# impala 0.2.0
* add `calibHier()` for hierarchical calibration, where each `theta_ind` group
  gets its own parameter drawn from N(theta0, Sigma0), with hyperpriors set by
  the new `setHierPriors()` (ported from python impala)
* add `calibClust()` for clustered calibration, where a Dirichlet process
  groups the `theta_ind` parameters into shared clusters, configured with the
  new `setClusterPriors()` (ported from python impala)
* add observation weights: `addVecExperiments(wt = )` raises each
  observation's likelihood contribution to the power of its weight. Implemented
  for `ModelBassPca_func()`; custom models opt in with `supports_wt = TRUE`.
  The Gibbs error variance update uses the weights exactly
* `theta_ind` in `addVecExperiments()` is now 1-based (default `rep(1, n)`) and
  validated, matching `s2_ind`
* emulator `lik_cov_inv()` methods accept `inds` to build the covariance of a
  subset of the observations
* `calibPool()` now honors `s2 = 'fix'`, holding the error variance at
  `sd_est^2` instead of sampling it by Metropolis-Hastings
* add `fixTheta()`/`unfixTheta()` to hold calibration parameters at fixed values,
  supporting conditional and cut-Bayes inference
* document that `sd_est`, `s2_df` and `s2_ind` in `addVecExperiments()` are
  indexed by measurement error group, so a separate standard deviation can be
  learned per group -- including one per component of `yobs` via
  `s2_ind = seq_along(yobs)`
* `addVecExperiments()` now validates `sd_est`, `s2_df` and `s2_ind` instead of
  silently dropping observations or corrupting the error prior; 0-based `s2_ind`
  (as used by python impala) is rejected with a message rather than leaving an
  empty group
* add optional `sd_lower`/`sd_upper` bounds on the learned measurement error
  standard deviations, enforced in both the Gibbs and Metropolis-Hastings `s2`
  updates
* fix `calibPool()` failing with "incorrect number of dimensions" for emulator
  models at a single temperature
* add ability to specify priors for theta (#8)
* cleaned up documentation
* numerical bugfixes
* added test harness (#11)

# impala 0.1.4
* Bugfixes for tempering swaps

# impala 0.1.3
* Bugfixes for dimensions

# impala 0.1.2
* Bugfixes to discrepancy

# impala 0.1.1
* Initial CRAN submission.
