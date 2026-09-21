# impala (development version)
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
