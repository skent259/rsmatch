
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsmatch

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/skent259/rsmatch.svg?branch=master)](https://travis-ci.com/skent259/rsmatch)
<!-- badges: end -->

The R package rsmatch implements various forms of **matching on
time-varying observational studies** to help in causal inference. The
package utilized the idea of matching through a series of risk sets
consisting of all subjects who are still “at risk” of treatment.
Subjects in this set treated at the next time point can be matched to
controls who are not treated to avoid bias in estimating treatment
effects.

Currently, we have methods for:

- **Balanced Risk Set Matching** with `brsmatch()`, based on the work of
  Li, Propert, and Rosenbaum (2001)
- **Propensity Score Matching with Time-Dependent Covariates** with
  `coxph_match()`, based on the work of Lu (2005)

## Installation

You can install the released version of rsmatch from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rsmatch")
```

Alternatively, you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("skent259/rsmatch")
```

## Example usage

Consider the `oasis` data which consists of 51 patients from the [Open
Access Series of Imaging Studies](https://www.oasis-brains.org/) that
are classified at each time point as having Alzheimer’s disease (AD) or
not. We can match subjects who are similar in their gender, education
level, socioeconomic status, age, and other MRI measurements, but one of
which moved from a cognitively normal state to AD in the next time
period and the other who remained cognitively normal.

``` r
library(rsmatch)
data(oasis)

pairs <- brsmatch(
  n_pairs = 5,
  data = oasis,
  id = "subject_id", time = "visit", trt_time = "time_of_ad",
  balance = TRUE, balance_covariates = c("m_f", "age")
)

na.omit(pairs)
#>    subject_id pair_id type
#> 5   OAS2_0014       1  trt
#> 15  OAS2_0043       2  all
#> 22  OAS2_0079       2  trt
#> 32  OAS2_0112       3  all
#> 36  OAS2_0124       1  all
#> 40  OAS2_0140       5  all
#> 41  OAS2_0150       3  trt
#> 43  OAS2_0160       4  trt
#> 49  OAS2_0182       4  all
#> 50  OAS2_0184       5  trt
```

These 10 subjects form the 5 best pairs according to balanced risk set
matching, where we consider pairs close if they have similar covariates
and we try to find balance on gender and age. We can see that the first
pair match closely on their covariates:

``` r
first_pair <- pairs[which(pairs$pair_id == 1), "subject_id"]
oasis[which(oasis$subject_id %in% first_pair), ]
#>    subject_id visit time_of_ad m_f educ ses age mr_delay e_tiv n_wbv   asf
#> 11  OAS2_0014     1          2   M   16   3  76        0  1602 0.697 1.096
#> 12  OAS2_0014     2          2   M   16   3  77      504  1590 0.696 1.104
#> 80  OAS2_0124     1         NA   M   16   3  70        0  1463 0.749 1.200
#> 81  OAS2_0124     2         NA   M   16   3  71      472  1479 0.750 1.187
```

We can also find matches using Propensity Score Matching with
Time-Dependent Covariates by Lu (2005)

``` r
pairs <- coxpsmatch(
  n_pairs = 5,
  data = oasis,
  id = "subject_id", time = "visit", trt_time = "time_of_ad"
)

first_pair <- pairs[which(pairs$pair_id == 1), "subject_id"]
oasis[which(oasis$subject_id %in% first_pair), ]
#> # A tibble: 4 × 11
#>   subject_id visit time_of_ad m_f    educ ses     age mr_delay e_tiv n_wbv   asf
#>   <chr>      <int>      <dbl> <chr> <int> <fct> <int>    <int> <int> <dbl> <dbl>
#> 1 OAS2_0009      1         NA M        12 2        68        0  1457 0.806  1.20
#> 2 OAS2_0009      2         NA M        12 2        69      576  1480 0.791  1.19
#> 3 OAS2_0046      1          2 F        15 2        83        0  1476 0.75   1.19
#> 4 OAS2_0046      2          2 F        15 2        85      575  1483 0.748  1.18
```

This picks a different first pair, but they are also close in
covariates.

## Citation

To cite package ‘rsmatch’ in publications use:

Kent S, Paukner M (2023). *rsmatch: Matching Methods for Time-varying
Observational Studies*. R package version 0.1.0,
<https://github.com/skent259/rsmatch>.

A BibTeX entry for LaTeX users is

``` bibtex
@Manual{,
  title = {rsmatch: Matching Methods for Time-varying Observational Studies},
  author = {Sean Kent and Mitchell Paukner},
  year = {2023},
  note = {R package version 0.1.0},
  url = {https://github.com/skent259/rsmatch},
}
```

## References

Li, Yunfei Paul, Kathleen J Propert, and Paul R Rosenbaum. 2001.
“Balanced Risk Set Matching.” Journal of the American Statistical
Association 96 (455): 870–82.

Lu, Bo. 2005. “Propensity Score Matching with Time-Dependent
Covariates.” Biometrics 61 (3): 721–28.
