# rsmatch 0.2.1

* `brsmatch()` and `coxpsmatch()` now handle NA rows via removing them (#16)
* Fix `brsmatch()` to iteratively reduce the number of pairs if too many are specified (#11)

# rsmatch 0.2.0

* Initial CRAN submission.
* BREAKING: rename `coxph_match()` to `coxpsmatch()`
  * Also changes to input arguments
  * Add exact matching functionality
* BREAKING: major changes to `brsmatch()`, including
  * Changes to input arguments
  * Add exact matching functionality
* Add `oasis` data for simple examples of matching
* Add README.md file
* Other documentation and internal improvements

# rsmatch 0.1.0

* Initial package release, including two main functions
* `brsmatch()` which implements Balanced Risk Set Matching based on the work of Li et al. (2001)
* `coxpsmatch()` which implements Propensity Score Matching with Time-Dependent Covariates based on the work of Lu (2005)
