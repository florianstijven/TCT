
# TCT

The goal the TCT package is to provide a set of easy to use functions to
perform time components test.

## Installation You can install the development

version of TCT from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("florianstijven/TCT")
```

## Example

This is a basic example which shows you how to use the Meta-TCT
methodology in combination with an estimated MMRM model. First the MMRM
is estimated. All we need from this estimated model is the estimated
mean vector and estimated variance-covariance matrix for this estimated
vector.

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
```

After fitting a MMRM model, we can use the `TCT_meta()` function to
perform the time-component test. The first step is to estimate the
slowing factor for each time point. Different approaches for inference
are available. Note that the estimators in these different approaches
are equivalent, but the corresponding measures of uncertainty are not
equivalent. Information on these approaches can be found in the function
documentation.

``` r
# TCT_meta() is the main function for performing time component tests. The estimated
# mean parameters of the MMRM model and associated variance-covariance matrix
# form the basis for the time component test.
TCT_fit = TCT_meta(
  time_points = 0:4,
  ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
  exp_estimates = coef(mmrm_fit)[5:8],
  vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
  interpolation = "spline",
  inference = "score",
  B = 1e4
)
# Summary of the results of the time component test. Here we have estimated a
# treatment effect parameter on the time scale for every measurement occasion.
summary(TCT_fit)
```

    ## Meta-Time Component Test: 
    ## 
    ## Coefficients: 
    ##               Value Std. Error  z value  p value                  CI
    ## arm_time1:2 0.73275   0.620432 -0.59571 0.551368 (-0.28696, 1.35407)
    ## arm_time1:3 0.87347   0.101620 -1.23894 0.215366 ( 0.66403, 1.07477)
    ## arm_time1:4 0.79370   0.080407 -2.13016 0.033159 ( 0.64368, 0.97891)
    ## arm_time1:5 0.75102   0.096581 -1.71943 0.085537 ( 0.60176, 1.03671)
    ##                  CI (bootstrap)
    ## arm_time1:2 (-0.27639, 1.35274)
    ## arm_time1:3 ( 0.66452, 1.07347)
    ## arm_time1:4 ( 0.64222, 0.97614)
    ## arm_time1:5 ( 0.60127, 1.03710)
    ## alpha = 0.05
    ## 
    ## Interpolation Method: spline
    ## Estimation and Inference: score

The next step is to estimate a common slowing factor. As before,
different approaches to inference are available. However, the estimators
in these different approaches are no longer equivalent. In the next code
chunk, score-based weighted with adaptively selected weights is
illustrated. These weight the different measurement occasions such the
the variance of the corresponding estimator is minimized.

``` r
# The object returned by TCT() can be used as argument to the TCT_common() function.
# This function estimates a single common slowing factor.
TCT_common_fit = TCT_meta_common(TCT_fit,
                                 B = 0,
                                 inference = "score",
                                 type = "custom")
summary(TCT_common_fit)
```

    ## Meta-Time Component Test - Common Acceleration Factor:
    ## 
    ## Estimated Common Acceleration Factor: 
    ##   Estimate Std. Error z value  p value                 CI
    ## z   0.7937   0.080407 -2.1302 0.033159 (0.64368, 0.97891)
    ## alpha = 0.05
    ## 
    ## Interpolation Method: spline
    ## Time Points Used in Estimator: 1 2 3 4
    ## Estimation and Inference: score
    ##  Type of Score Test: custom
    ##  Weights: 0 0 1 0
    ## 
    ## Test for proportional slowing factor:
    ##   Df  Chisq p.value
    ## 1  3 4.1219  0.2486
