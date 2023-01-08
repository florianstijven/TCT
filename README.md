
# TCT

<!-- badges: start -->
<!-- badges: end -->

The goal the TCT package is to provide a set of easy to use functions to perform 
time components test. 

## Installation

You can install the development version of TCT from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("florianstijven/TCT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TCT)
library(dplyr)
# Example data set transformed to format required by TCT()
data_test = simulated_test_trial %>%
  mutate(
    time_int = (Week %/% 25) + 1,
    arm_time = ifelse(time_int == 1L,
                      "baseline",
                      paste0(arm, ":", time_int))
  )
# Fit a MMRM model to the data set. The parameter estimates of this model form 
# the basis to perform the time component tests.
mmrm_fit = analyze_mmrm(data_test)

# TCT() is the main function for performing time component tests. The estimated
# mean parameters of the MMRM model and associated variance-covariance matrix
# form the basis for the time component test.
TCT_fit = TCT(time_points = 0:4,
              ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)], 
              exp_estimates = coef(mmrm_fit)[5:8], 
              vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)], 
              interpolation = "spline", 
              B = 0)
# Summary of the results of the time component test.
summary(x = TCT_fit)
# If one assumes proportional 
TCT_common_fit = TCT_common(TCT_fit, B = 1e4, TRUE)
summary(TCT_common_fit)

TCT_fit = TCT(time_points = 0:4,
              ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)], 
              exp_estimates = coef(mmrm_fit)[5:8], 
              vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)], 
              interpolation = "linear", 
              B = 0)
summary(x = TCT_fit)
TCT_common_fit = TCT_common(TCT_fit, B = 1e4, bs_fix_vcov = TRUE)
summary(TCT_common_fit)

```

