test_that("Optimisation for weights vector works", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")
  # t-value for TCT score test
  weights_opt1 = optimize_weights(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    ref_fun = ref_fun,
    interpolation = "spline",
    vcov = vcov_mmrm,
    j = 1:4
  )
  weights_opt2 = optimize_weights(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    ref_fun = ref_fun,
    interpolation = "spline",
    vcov = diag(1, nrow = 9),
    j = 1:4
  )
  output_vector = c(weights_opt1, weights_opt2)
  check_vector = c(
    0,
    0,
    1,
    0,
    0.025362298,
    0.240318869,
    0.387036732,
    0.347282101
  )
  expect_equal(output_vector,
               check_vector)
})
