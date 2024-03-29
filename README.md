
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/imbi-heidelberg/blindrecalc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/imbi-heidelberg/blindrecalc/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/imbi-heidelberg/blindrecalc/branch/master/graph/badge.svg)](https://app.codecov.io/gh/imbi-heidelberg/blindrecalc?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/blindrecalc)](https://cran.r-project.org/package=blindrecalc)
<!-- badges: end -->

# blindrecalc

blindrecalc facilitates the planning of a clinical trial with an
internal pilot study and blinded sample size recalculation.

## Installation

Install the current CRAN version of blindrecalc with:

``` r
install.packages("blindrecalc")
```

Or install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("imbi-heidelberg/blindrecalc")
```

## Usage

blindrecalc currently supports continuous and binary endpoints for
superiority and non-inferiority test problems. Continuous endpoints are
analyzed using Student’s t-test, binary endpoints are analyzed using the
Chi-squared test for superiority trials and the Farrington-Manning test
for non-inferiority trials. Each design can be defined using a
setup-function: `setupStudent`, `setupChiSquare` and
`setupFarringtonManning`. For example, to setup a superiority trial with
a continuous endpoint:

``` r
library(blindrecalc)
design <- setupStudent(alpha = 0.025, beta = 0.2, r = 1, delta = 5)
```

`alpha` and `beta` refer to the type 1 and type 2 error rate, `r` is the
sample size allocation ratio and `delta`is the effect size between the
null and the alternative hypothesis. For a non-inferiority trial with a
shifted t-test, additionally the argument `delta_NI` must be specified.

To calculate the sample size for a fixed design, use `n_fix`:

``` r
n_fix(design, nuisance = c(5, 10, 15))
#> [1]  31.39552 125.58208 282.55967
```

`nuisance` refers to the nuisance parameter of the design, which in the
case of the t-test is the common variance of the outcome variable.

To calculate the type 1 error rate of the design using blinded sample
size recalculation, use `toer`:

``` r
toer(design, n1 = c(30, 60, 90), nuisance = 10, recalculation = TRUE)
#> [1] 0.0259 0.0235 0.0252
```

`n1` refers to the sample size of the internal pilot study
`recalculation = TRUE` specifices that the type 1 error rate for a
design with blinded sample size recalculation should be computed.

To compute the power of the design, use `pow`:

``` r
pow(design, n1 = c(30, 60, 90), nuisance = 10, recalculation = TRUE)
#> [1] 0.7877 0.8039 0.8056
```

To calculate the distribution of the total sample sizes use `n_dist`:

``` r
n_dist(design, n1 = c(30, 60, 90), nuisance = 10)
#>     n_1 = 30        n_1 = 60      n_1 = 90    
#>  Min.   : 36.0   Min.   : 63   Min.   : 90.0  
#>  1st Qu.:109.0   1st Qu.:117   1st Qu.:120.0  
#>  Median :131.0   Median :132   Median :133.0  
#>  Mean   :134.1   Mean   :134   Mean   :134.5  
#>  3rd Qu.:155.0   3rd Qu.:150   3rd Qu.:147.0  
#>  Max.   :322.0   Max.   :282   Max.   :223.0
```

## Reference

A paper describing blindrecalc can be found
[here](https://journal.r-project.org/articles/RJ-2022-001/).
