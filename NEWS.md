# impala (development version)
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
