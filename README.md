# rsmatch
<!-- badges: start -->
[![Travis build status](https://travis-ci.com/skent259/rsmatch.svg?branch=master)](https://travis-ci.com/skent259/rsmatch)
<!-- badges: end -->

This repo is an R package that implements various forms of **matching on time-varying observational studies** to help in causal inference.  
The work is in progress and we hope to expand the documentation shortly.

Currently, we have methods for:

- **Balanced Risk Set Matching** with `brsmatch()`, based on the work of Li et al. (2001) 
- **Propensity Score Matching with Time-Dependent Covariates** with `coxph_match()`, based on the work of Lu (2005)




## References 

Li, Yunfei Paul, Kathleen J Propert, and Paul R Rosenbaum. 2001. “Balanced Risk Set Matching.” Journal of the American Statistical Association 96 (455): 870–82. https://doi.org/10.1198/016214501753208573.

Lu, Bo. 2005. “Propensity Score Matching with Time-Dependent Covariates.” Biometrics 61 (3): 721–28.
