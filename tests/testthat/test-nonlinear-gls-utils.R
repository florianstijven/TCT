test_that("numerical and analytical derivatives match for GLS criterion function", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
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
  # Run the nonlinear GLS estimator.
  criterion_function = nonlinear_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline"
  )
  gradient_function = gradient_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline"
  )
  x = c(22.3683499924, 23.8805436692, 26.2779144133, 29.4903392124, 31.6484238040, 0.7758874894)
  output_vctr = as.numeric(gradient_function(x))
  check_vctr = numDeriv::grad(func = criterion_function,
                              x = x)
  expect_equal(output_vctr, check_vctr)
})

test_that("numerical and analytical derivatives match for GLS criterion function (reduced model)", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
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
  # Run the nonlinear GLS estimator.
  criterion_function = nonlinear_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    gamma_0 = 1.1
  )
  gradient_function = gradient_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    gamma_0 = 1.1
  )
  x = c(22.3683499924, 23.8805436692, 26.2779144133, 29.4903392124, 31.6484238040)
  output_vctr = as.numeric(gradient_function(x))
  check_vctr = numDeriv::grad(func = criterion_function,
                              x = x)
  expect_equal(output_vctr, check_vctr)
})
