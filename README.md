
# TCT The goal the TCT package is to provide a set of easy to use functions to

perform time components test. \## Installation You can install the
development version of TCT from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("florianstijven/TCT")
```

## Example This is a basic example which shows you how to use the time-component

test methodology in combination with an estimated MMRM model. First the
MMRM is estimated. All we need from this estimated model is the
estimated mean vector and estimated variance-covariance matrix for this
estimated vector.

``` r
library(TCT)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
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

After fitting a MMRM model, we can use the `TCT()` function to perform
the time-component test. The first step is to estimate the slowing
factor for each time point.

``` r
# TCT() is the main function for performing time component tests. The estimated
# mean parameters of the MMRM model and associated variance-covariance matrix
# form the basis for the time component test.
TCT_fit = TCT(time_points = 0:4,
              ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
              exp_estimates = coef(mmrm_fit)[5:8],
              vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
              interpolation = "spline",
              B = 1e4)
# Summary of the results of the time component test. Here we have estimated a
# treatment effect summary(x = TCT_fit) parameter on the time scale for every
# measurement occasion. 
summary(TCT_fit)
```

    ## Time Component Test: time-based treatment effects
    ## 
    ## Coefficients: 
    ##               Value Std. Error (delta) z-value (delta) p-value (delta)
    ## arm_time1:2 0.73268           0.620592         0.43075       0.6666514
    ## arm_time1:3 0.87351           0.101608         1.24491       0.2131653
    ## arm_time1:4 0.79369           0.080397         2.56612       0.0102844
    ## arm_time1:5 0.75101           0.096568         2.57841       0.0099256
    ##                      CI (delta)     CI (bootstrap)
    ## arm_time1:2 (-0.48366, 1.94902) (-0.28498, 1.3522)
    ## arm_time1:3 ( 0.67436, 1.07266) ( 0.66815, 1.0731)
    ## arm_time1:4 ( 0.63612, 0.95127) ( 0.64753, 0.9755)
    ## arm_time1:5 ( 0.56174, 0.94028) ( 0.60534, 1.0356)
    ## alpha = 0.05
    ##  Interpolation Method: spline

The next step is to estimate a common slowing factor.

``` r
# The object returned by TCT() can be used as argument to the TCT_common() function.
# This function estimates a single common slowing factor.
TCT_common_fit = TCT_common(TCT_fit, B = 1e4, TRUE)
summary(TCT_common_fit)
```

    ## Time Component Test: proportional slowing
    ## 
    ## Coefficients: 
    ##     Value Std. Error (delta) z-value (delta) p-value (delta)         CI (delta)
    ## 1 0.80275           0.067897          2.9051       0.0036714 (0.66968, 0.93583)
    ##      CI (bootstrap) p-value (bootstrap)
    ## 1 (0.68494, 0.9832)               0.035
    ## alpha = 0.05
    ## 
    ## Test for proportional slowing factor:
    ##   Df  Chisq p.value
    ## 1  3 4.1253 0.24825
