test_that("numerical and analytical derivatives match for GLS criterion function", {
  # Run the nonlinear GLS estimator.
  criterion_function = nonlinear_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  gradient_function = gradient_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  x = c(22.3683499924, 23.8805436692, 26.2779144133, 29.4903392124, 31.6484238040, 0.7758874894)
  output_vctr = as.numeric(gradient_function(x))
  check_vctr = numDeriv::grad(func = criterion_function,
                              x = x)
  expect_equal(output_vctr, check_vctr)
})

test_that("numerical and analytical derivatives match for GLS criterion function (reduced model)", {
  # Run the nonlinear GLS estimator.
  criterion_function = nonlinear_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    gamma_0 = 1.1
  )
  gradient_function = gradient_gls_criterion_constructor(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    gamma_0 = 1.1
  )
  x = c(22.3683499924, 23.8805436692, 26.2779144133, 29.4903392124, 31.6484238040)
  output_vctr = as.numeric(gradient_function(x))
  check_vctr = numDeriv::grad(func = criterion_function,
                              x = x)
  expect_equal(output_vctr, check_vctr)
})

