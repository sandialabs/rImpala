## Resubmission
This is a resubmission of impala 0.2.0. The first submission passed the
incoming pre-tests (Status: OK on Windows and Debian) but the Windows check
took over 10 minutes, mostly in the test suite (453s). The long-chain MCMC
posterior-accuracy tests now use `testthat::skip_on_cran()`; they still run
locally and on CI. Tests now take about 12s locally (previously 77s).

## Test environments
* windows, os x, ubuntu 24 (on githubactions), R 4.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is an update to an existing CRAN package (0.1.4 -> 0.2.0).
